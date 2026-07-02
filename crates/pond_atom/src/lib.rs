//! # pond_atom — Per-point peridynamics data, materials, and point clouds
//!
//! This crate provides the per-point data structures, material system, and
//! reference-configuration point-cloud generation for peridynamics (PD)
//! simulations in POND:
//!
//! - **[`PdAtom`]** — Per-point extension data (reference position, point
//!   volume, scalar damage) registered as a SOIL [`AtomData`] column. The
//!   point `volume` is `#[forward]` so ghost points carry it across MPI ranks
//!   (a bond force on point `i` needs the *partner's* volume `Vⱼ`).
//! - **[`PdMaterialTable`]** — Named PD materials. Each material stores the
//!   bulk modulus `K`, the constant bond micromodulus `c = 18K/(πδ⁴)`, and the
//!   critical bond stretch `s₀`. Bond-based PD fixes Poisson's ratio at ν = ¼,
//!   so `K = 2E/3`.
//! - **[`PdLatticePlugin`]** — Generates a regular point cloud from
//!   `[[pond.lattice]]` config blocks: a uniform grid of points, each carrying
//!   a cell volume `Δx³`, mass `ρΔx³`, and a neighbour cutoff equal to the
//!   horizon δ so SOIL's neighbour machinery sees the whole family.
//!
//! ## TOML configuration
//!
//! ```toml
//! [pond]
//! horizon = 0.0015          # PD horizon δ [m] (typically ≈ 3·spacing)
//!
//! [[pond.materials]]
//! name = "pmma"
//! youngs_mod = 3.0e9        # Young's modulus E [Pa]
//! density = 1190.0          # [kg/m³]
//! critical_stretch = 0.04   # s₀ (omit to derive from fracture_energy, or both to stay elastic)
//! # fracture_energy = 50.0  # G₀ [J/m²]  → s₀ = sqrt(5 G₀ / (9 K δ))
//!
//! [[pond.lattice]]
//! material = "pmma"
//! min = [0.0, 0.0, 0.0]     # reference-box lower corner [m]
//! max = [0.05, 0.01, 0.001] # reference-box upper corner [m]
//! spacing = 0.0005          # point spacing Δx [m]
//! # velocity = [0.0, 0.0, 0.0]  # optional initial velocity [m/s]
//! ```

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use serde::Deserialize;
use soil_derive::AtomData; // derive macro (the trait of the same name comes from soil_core)

use soil_core::{
    register_atom_data, Atom, AtomData, AtomPlugin, CommResource, Config, Domain, ScheduleSetupSet,
};

/// Poisson's ratio is fixed at ¼ in the bond-based peridynamic model.
pub const BOND_BASED_POISSON: f64 = 0.25;

// ── Config structs ───────────────────────────────────────────────────────────

/// The `[pond]` configuration section: horizon plus material and lattice lists.
#[derive(Deserialize, Clone, Default)]
pub struct PondConfig {
    /// Peridynamic horizon δ [m] — the bond interaction radius. A point bonds to
    /// every other point within δ in the reference configuration.
    pub horizon: f64,
    /// Material definitions (`[[pond.materials]]`).
    pub materials: Option<Vec<PdMaterialDef>>,
    /// Reference point-cloud blocks (`[[pond.lattice]]`).
    pub lattice: Option<Vec<LatticeDef>>,
}

/// A single `[[pond.materials]]` entry.
#[derive(Deserialize, Clone)]
pub struct PdMaterialDef {
    /// Material name, referenced by `[[pond.lattice]] material`.
    pub name: String,
    /// Young's modulus E [Pa].
    pub youngs_mod: f64,
    /// Mass density ρ [kg/m³].
    pub density: f64,
    /// Critical bond stretch s₀ (dimensionless). Takes precedence over
    /// `fracture_energy` when both are given.
    pub critical_stretch: Option<f64>,
    /// Fracture energy G₀ [J/m²]. Converted to a critical stretch via
    /// `s₀ = sqrt(5 G₀ / (9 K δ))` (3-D constant-micromodulus model).
    pub fracture_energy: Option<f64>,
}

/// A single `[[pond.lattice]]` point-cloud block.
#[derive(Deserialize, Clone)]
pub struct LatticeDef {
    /// Name of the material these points are made of.
    pub material: String,
    /// Reference-box lower corner [m].
    pub min: [f64; 3],
    /// Reference-box upper corner [m].
    pub max: [f64; 3],
    /// Point spacing Δx [m]. The cell volume per point is `Δx³`.
    pub spacing: f64,
    /// Optional uniform initial velocity [m/s].
    pub velocity: Option<[f64; 3]>,
}

// ── Material table ────────────────────────────────────────────────────────────

/// A resolved PD material: density plus the precomputed bond constants for the
/// configured horizon.
#[derive(Clone, Debug)]
pub struct PdMaterial {
    /// Material name.
    pub name: String,
    /// Mass density ρ [kg/m³].
    pub density: f64,
    /// Young's modulus E [Pa].
    pub youngs_mod: f64,
    /// Bulk modulus K = E / (3(1 − 2ν)) with ν = ¼, i.e. K = 2E/3 [Pa].
    pub bulk_modulus: f64,
    /// Constant bond micromodulus c = 18K / (πδ⁴) [N/m⁶].
    pub micromodulus: f64,
    /// Critical bond stretch s₀ (dimensionless). `INFINITY` ⇒ unbreakable
    /// (purely elastic).
    pub critical_stretch: f64,
}

/// Named PD materials with their precomputed bond constants. Added as a
/// resource by [`PdAtomPlugin`].
#[derive(Clone, Debug, Default)]
pub struct PdMaterialTable {
    /// The horizon δ [m] these constants were built for.
    pub horizon: f64,
    /// Resolved materials, indexed by atom type.
    pub materials: Vec<PdMaterial>,
}

impl PdMaterialTable {
    /// Builds the table from the parsed `[pond]` config, precomputing the
    /// bond constants for the given horizon.
    pub fn from_config(cfg: &PondConfig) -> Self {
        let horizon = cfg.horizon;
        assert!(horizon > 0.0, "pond.horizon must be > 0");

        let mut materials = Vec::new();
        if let Some(defs) = &cfg.materials {
            for def in defs {
                // Bond-based PD fixes ν = ¼  ⇒  K = E / (3(1 − 2ν)) = 2E/3.
                let k = def.youngs_mod / (3.0 * (1.0 - 2.0 * BOND_BASED_POISSON));
                let c = 18.0 * k / (std::f64::consts::PI * horizon.powi(4));
                let s0 = match (def.critical_stretch, def.fracture_energy) {
                    (Some(s0), _) => s0,
                    (None, Some(g0)) => (5.0 * g0 / (9.0 * k * horizon)).sqrt(),
                    (None, None) => f64::INFINITY, // elastic, never breaks
                };
                materials.push(PdMaterial {
                    name: def.name.clone(),
                    density: def.density,
                    youngs_mod: def.youngs_mod,
                    bulk_modulus: k,
                    micromodulus: c,
                    critical_stretch: s0,
                });
            }
        }
        PdMaterialTable { horizon, materials }
    }

    /// Index of a material by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.materials.iter().position(|m| m.name == name)
    }

    /// Effective micromodulus for a bond between two material types — the
    /// arithmetic mean of the two single-material micromoduli.
    #[inline]
    pub fn micromodulus_ij(&self, ti: usize, tj: usize) -> f64 {
        0.5 * (self.materials[ti].micromodulus + self.materials[tj].micromodulus)
    }

    /// Critical stretch for a bond between two material types — the smaller of
    /// the two (the bond fails at the weaker endpoint).
    #[inline]
    pub fn critical_stretch_ij(&self, ti: usize, tj: usize) -> f64 {
        self.materials[ti]
            .critical_stretch
            .min(self.materials[tj].critical_stretch)
    }
}

// ── Per-point data ────────────────────────────────────────────────────────────

/// Per-point peridynamic extension data, registered as a SOIL [`AtomData`]
/// column.
///
/// - `ref_pos` — the **reference** (undeformed) position. Reference *bond
///   lengths* live in the [`soil_core::BondStore`] entries, so `ref_pos` is not
///   communicated; it is kept locally for boundary conditions and diagnostics.
/// - `volume` — the cell volume Δx³ represented by the point. Marked
///   `#[forward]` so a ghost point carries its volume to the rank that needs it
///   when summing `Vⱼ` into a bond force.
/// - `damage` — the fraction of the point's bonds that have broken, in `[0, 1]`.
///   This is the field a crack shows up in.
#[derive(AtomData)]
pub struct PdAtom {
    /// Reference (undeformed) position [m].
    pub ref_pos: Vec<[f64; 3]>,
    /// Cell volume represented by the point [m³].
    #[forward]
    pub volume: Vec<f64>,
    /// Initial bond count of the point (the damage denominator). Set once in
    /// the reference configuration when the family is built.
    pub n0: Vec<f64>,
    /// Local damage: fraction of broken bonds, in `[0, 1]`.
    pub damage: Vec<f64>,
}

impl Default for PdAtom {
    fn default() -> Self {
        Self::new()
    }
}

impl PdAtom {
    /// Creates an empty `PdAtom` with no points. Columns are appended during
    /// lattice generation.
    pub fn new() -> Self {
        PdAtom {
            ref_pos: Vec::new(),
            volume: Vec::new(),
            n0: Vec::new(),
            damage: Vec::new(),
        }
    }
}

// ── Plugins ───────────────────────────────────────────────────────────────────

/// Registers [`PdAtom`] and builds the [`PdMaterialTable`] from `[pond]` config.
pub struct PdAtomPlugin;

impl Plugin for PdAtomPlugin {
    fn provides(&self) -> Vec<&str> {
        vec!["pd_points"]
    }

    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Peridynamic horizon and materials
[pond]
horizon = 0.0015

[[pond.materials]]
name = "pmma"
youngs_mod = 3.0e9
density = 1190.0
critical_stretch = 0.04"#,
        )
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(AtomPlugin);

        register_atom_data!(app, PdAtom::new());

        let cfg = Config::load::<PondConfig>(app, "pond");
        let table = PdMaterialTable::from_config(&cfg);

        app.add_resource(cfg);
        app.add_resource(table);
        app.add_setup_system(set_pd_ntypes, ScheduleSetupSet::Setup);
    }
}

/// Setup system that sets `Atom::ntypes` from the number of registered materials.
fn set_pd_ntypes(mut atoms: ResMut<Atom>, table: Res<PdMaterialTable>) {
    atoms.ntypes = table.materials.len().max(1);
}

/// Generates the reference point cloud from `[[pond.lattice]]` blocks.
pub struct PdLatticePlugin;

impl Plugin for PdLatticePlugin {
    fn requires(&self) -> Vec<&str> {
        vec!["pd_points"]
    }

    fn build(&self, app: &mut App) {
        // Must run after domain decomposition so each rank's subdomain bounds
        // are populated before points are placed and filtered to the subdomain.
        app.add_setup_system(
            pond_lattice_insert.after("domain_read_input"),
            ScheduleSetupSet::Setup,
        );
    }
}

/// Whether `pos` lies inside this rank's subdomain, using a half-open interval
/// `[sub_domain_low, sub_domain_high)` per axis.
#[inline]
fn owns_position(domain: &Domain, pos: &[f64; 3]) -> bool {
    (0..3).all(|d| pos[d] >= domain.sub_domain_low[d] && pos[d] < domain.sub_domain_high[d])
}

/// Setup system: lays down a uniform grid of points for every `[[pond.lattice]]`
/// block, placing each point at a cell centre and giving it a cell volume,
/// mass, horizon cutoff, and the configured material type.
fn pond_lattice_insert(
    mut atoms: ResMut<Atom>,
    registry: Res<soil_core::AtomDataRegistry>,
    cfg: Res<PondConfig>,
    table: Res<PdMaterialTable>,
    domain: Res<Domain>,
    comm: Res<CommResource>,
) {
    let lattices = match &cfg.lattice {
        Some(l) if !l.is_empty() => l,
        _ => return,
    };

    let horizon = table.horizon;
    let mut pd = registry.expect_mut::<PdAtom>("pond_lattice_insert");

    // Tags are assigned globally and sequentially. For the single-process slice
    // this is a simple running counter; cross-rank tag uniqueness is a TODO for
    // the parallel path.
    let mut next_tag: u32 = atoms.natoms as u32 + 1;
    let mut placed: u64 = 0;

    for lat in lattices {
        let mat_idx = table.index_of(&lat.material).unwrap_or_else(|| {
            panic!("pond.lattice references unknown material '{}'", lat.material)
        });
        let density = table.materials[mat_idx].density;
        let dx = lat.spacing;
        assert!(dx > 0.0, "pond.lattice spacing must be > 0");

        let volume = dx * dx * dx;
        let mass = density * volume;
        let inv_mass = 1.0 / mass;
        let vel = lat.velocity.unwrap_or([0.0; 3]);

        // Number of cells along each axis; points sit at cell centres so the
        // reference box [min, max] is tiled by Δx³ cells.
        let n: [i64; 3] = std::array::from_fn(|d| {
            (((lat.max[d] - lat.min[d]) / dx).round() as i64).max(0)
        });

        for ix in 0..n[0] {
            for iy in 0..n[1] {
                for iz in 0..n[2] {
                    let pos = [
                        lat.min[0] + (ix as f64 + 0.5) * dx,
                        lat.min[1] + (iy as f64 + 0.5) * dx,
                        lat.min[2] + (iz as f64 + 0.5) * dx,
                    ];
                    let tag = next_tag;
                    next_tag += 1;

                    if !owns_position(&domain, &pos) {
                        continue;
                    }

                    // Base SOIL atom columns.
                    atoms.natoms += 1;
                    atoms.nlocal += 1;
                    atoms.tag.push(tag);
                    atoms.origin_index.push(0);
                    atoms.cutoff_radius.push(horizon);
                    atoms.image.push([0, 0, 0]);
                    atoms.is_ghost.push(false);
                    atoms.pos.push(pos);
                    atoms.vel.push(vel);
                    atoms.force.push([0.0; 3]);
                    atoms.mass.push(mass);
                    atoms.inv_mass.push(inv_mass);
                    atoms.atom_type.push(mat_idx as u32);

                    // PD columns. `n0` is filled when the family is built.
                    pd.ref_pos.push(pos);
                    pd.volume.push(volume);
                    pd.n0.push(0.0);
                    pd.damage.push(0.0);

                    placed += 1;
                }
            }
        }
    }

    let total = comm.all_reduce_sum_f64(placed as f64) as u64;
    if comm.rank() == 0 {
        println!("PdLattice: placed {total} peridynamic points");
    }
}
