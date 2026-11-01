//! # dev_soil_peri — Peridynamic Off-lattice Nonlocal Deformation
//!
//! dev_soil_peri is a plugin-based simulation framework for **peridynamics (PD)** built
//! on the GRASS → SOIL stack. This umbrella crate re-exports every dev_soil_peri plugin
//! crate and provides convenience plugin groups so a complete simulation reads
//! in a few lines.
//!
//! # Quick start
//!
//! ```rust,ignore
//! use peri_core::prelude::*;
//!
//! fn main() {
//!     let mut app = App::new();
//!     app.add_plugins(CorePlugins)
//!        .add_plugins(BondBasedDefaultPlugins);
//!     app.start();
//! }
//! ```
//!
//! The simulation reads its parameters from a TOML config file passed as a
//! command-line argument (e.g. `cargo run -- config.toml`).
//!
//! # Plugin groups
//!
//! | Group | Purpose |
//! |---|---|
//! | [`CorePlugins`] | Infrastructure: input, communication, domain, neighbour lists, groups, run loop, output |
//! | [`BondBasedDefaultPlugins`] | Bond-based PD: point cloud, velocity-Verlet integration, nonlocal bond force with fracture |
//!
//! # Crate organisation
//!
//! | Crate | Description |
//! |---|---|
//! | [`peri_atom`] | per-point PD data (`PeriAtom`), PD materials, lattice point-cloud generation |
//! | [`peri_bond`] | bond-based PD: constant-micromodulus force, critical-stretch breakage, damage |
//! | [`soil_core`] | substrate: atoms, domain decomposition, communication, neighbour lists, bonds |
//! | [`soil_verlet`] | velocity-Verlet integration |
//! | [`soil_print`] | thermo output, dump files, restart |

// --- PD crates ---

/// Per-point PD data, materials, and reference point-cloud generation.
pub use peri_atom;

/// Bond-based peridynamics: nonlocal pairwise force, fracture, damage.
pub use peri_bond;

// --- Infrastructure crates ---

/// Core simulation types: atoms, config, domain decomposition, communication, bonds.
pub use soil_core;

/// Velocity-Verlet time integration.
pub use soil_verlet;

/// Thermo output, dump files, and restart file I/O.
pub use soil_print;

use grass_app::prelude::*;

/// Core simulation infrastructure plugin group — SOIL's shared
/// [`soil_plugins::CorePlugins`], the same group DIRT and dev_soil_sph use.
/// Integration is **not** included; [`BondBasedDefaultPlugins`] adds
/// velocity-Verlet, which also applies `[[run]] dt`.
pub use soil_plugins::CorePlugins;

/// Bond-based peridynamics plugin group.
///
/// Includes, in registration order:
/// - [`PeriAtomPlugin`](peri_atom::PeriAtomPlugin) — `PeriAtom` data and `PeriMaterialTable`
/// - [`PeriLatticePlugin`](peri_atom::PeriLatticePlugin) — reference point cloud from `[[dev_soil_peri.lattice]]`
/// - [`VelocityVerletPlugin`](soil_verlet::VelocityVerletPlugin) — translational integration
/// - [`BondBasedPlugin`](peri_bond::BondBasedPlugin) — family construction + nonlocal bond force
pub struct BondBasedDefaultPlugins;

impl PluginGroup for BondBasedDefaultPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(peri_atom::PeriAtomPlugin)
            .add(peri_atom::PeriLatticePlugin)
            .add(soil_verlet::VelocityVerletPlugin::new())
            .add(peri_bond::BondBasedPlugin)
    }
}

/// The dev_soil_peri prelude — `use peri_core::prelude::*;` for everything a typical
/// simulation needs.
pub mod prelude {
    // Plugin groups defined in this crate.
    pub use crate::{BondBasedDefaultPlugins, CorePlugins};

    // PD plugins and types.
    pub use peri_atom::{
        LatticeDef, PeriAtom, PeriAtomPlugin, PeriLatticePlugin, PeriMaterial, PeriMaterialTable,
        PondConfig,
    };
    pub use peri_bond::{bond_based_force, BondBasedPlugin, BOND_BROKEN, BOND_INTACT};

    // Shared infrastructure.
    pub use soil_deform::{DeformConfig, DeformPlugin, DeformState};
    pub use soil_fixes::{GravityConfig, GravityPlugin, PinDef, PinRegistry, PinState, SoilFixesPlugin};

    // Derive macros.
    pub use grass_derive::{ScheduleSet, StageEnum};
    pub use soil_derive::AtomData;

    // Core framework re-exports (glob).
    pub use grass_app::prelude::*;
    pub use grass_scheduler::prelude::*;
    pub use soil_core::ParticleSimScheduleSet;
    pub use soil_core::*;
    pub use soil_print::*;
    pub use soil_verlet::*;
}
