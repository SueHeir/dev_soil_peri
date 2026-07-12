# peri_atom

Per-point peridynamics data, materials, and reference point-cloud generation
for [dev_soil_peri](../../README.md).

- **`PeriAtom`** — SOIL `AtomData` column carrying each point's reference position,
  cell volume (`#[forward]` so ghosts have it), and scalar damage.
- **`PeriMaterialTable`** — named PD materials with precomputed bond constants for
  the configured horizon: bulk modulus `K = 2E/3` (ν = ¼ in bond-based PD),
  micromodulus `c = 18K/(πδ⁴)`, and critical stretch `s₀` (given directly or
  derived from fracture energy `G₀`).
- **`PeriLatticePlugin`** — generates a uniform point cloud from `[[dev_soil_peri.lattice]]`
  blocks, one point per `Δx³` cell with mass `ρΔx³` and cutoff = horizon. Lattice
  tags are deterministic lattice-cell identities starting above the global
  pre-insertion maximum `u32` tag, so MPI ranks agree on halo partners while
  pre-existing atoms remain disjoint.

See the [workspace README](../../README.md) for the theory summary and TOML
configuration.
