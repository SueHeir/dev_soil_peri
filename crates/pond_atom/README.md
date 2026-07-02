# pond_atom

Per-point peridynamics data, materials, and reference point-cloud generation
for [POND](../../README.md).

- **`PdAtom`** — SOIL `AtomData` column carrying each point's reference position,
  cell volume (`#[forward]` so ghosts have it), and scalar damage.
- **`PdMaterialTable`** — named PD materials with precomputed bond constants for
  the configured horizon: bulk modulus `K = 2E/3` (ν = ¼ in bond-based PD),
  micromodulus `c = 18K/(πδ⁴)`, and critical stretch `s₀` (given directly or
  derived from fracture energy `G₀`).
- **`PdLatticePlugin`** — generates a uniform point cloud from `[[pond.lattice]]`
  blocks, one point per `Δx³` cell with mass `ρΔx³` and cutoff = horizon.

See the [workspace README](../../README.md) for the theory summary and TOML
configuration.
