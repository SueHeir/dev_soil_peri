# POND — Peridynamic Off-lattice Nonlocal Deformation

Peridynamics (PD) physics in the **GRASS → SOIL → {DIRT, POND}** stack.

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists
       ├─ DIRT   DEM physics: contact, bonds, walls, clumps
       └─ POND   peridynamics: nonlocal bond forces, fracture, damage   ← you are here
```

POND is the peridynamics sibling of [DIRT](https://github.com/SueHeir/dirt). Where
DIRT resolves short-range *contact* between discrete grains, POND resolves the
*nonlocal* interaction of continuum material points: every point interacts with
every neighbour inside a finite **horizon** δ through a bond, and material
failure emerges naturally as those bonds break. It rides the
[SOIL](https://github.com/SueHeir/soil) substrate for all the method-agnostic
machinery — atom data, domain decomposition, halo exchange, neighbour lists,
time integration — and adds only the peridynamic constitutive models.

The fit is close to exact: a PD point cloud *is* a SOIL atom set, the PD horizon
*is* a SOIL neighbour cutoff, and the PD family (the fixed set of bonds formed in
the reference configuration) *is* a SOIL [`BondStore`]. POND therefore reuses
SOIL's persistent, tag-addressed, MPI-aware bond storage rather than inventing
its own.

## Layout

`pond_core` is the batteries-included umbrella crate (depend on this):

```rust
use pond_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
       .add_plugins(BondBasedDefaultPlugins);
    app.start();
}
```

| crate | role |
|---|---|
| [`pond_core`](crates/pond_core/README.md) | umbrella: `CorePlugins`, `BondBasedDefaultPlugins`, prelude |
| [`pond_atom`](crates/pond_atom/README.md) | per-point PD data (`PdAtom`: reference position, volume, damage), PD materials, lattice point-cloud generation |
| [`pond_bond`](crates/pond_bond/README.md) | bond-based peridynamics: constant-micromodulus pairwise force, critical-stretch bond breakage, damage field |

## Theory in one paragraph

Bond-based peridynamics replaces the local stress divergence of classical
elasticity with an integral of pairwise bond forces over the horizon. For a
point **x** with reference bond **ξ = x′ − x** to a neighbour **x′** and relative
displacement **η**, the bond stretch is `s = (|ξ+η| − |ξ|) / |ξ|`, and the
pairwise force density is `f = c · s · (ξ+η)/|ξ+η|`. The micromodulus for the
3-D constant-`c` model is `c = 18K / (π δ⁴)`, with bulk modulus `K`. Bond-based
PD fixes Poisson's ratio at ν = 1/4 (hence `K = 2E/3`). A bond breaks
irreversibly once its stretch exceeds the critical stretch `s₀ = √(5G₀ / 9Kδ)`
(or a user-specified value); the local **damage** is the fraction of a point's
bonds that have broken, and a crack is simply the locus where damage → 1.

> **State-based PD** (ordinary and non-ordinary/correspondence), which removes
> the ν = 1/4 restriction, is the planned next tier (`pond_state`) on this same
> substrate.

## Building

Clone POND and build — its [soil](https://github.com/SueHeir/soil) and
[grass](https://github.com/SueHeir/grass) dependencies are pulled from GitHub
automatically.

```bash
git clone https://github.com/SueHeir/pond
cd pond
cargo run --release --example bb_bar_tension --no-default-features -- examples/bb_bar_tension/config.toml
```

The `--no-default-features` flag runs single-process (no MPI), which is the
supported configuration for the current bond-based slice — see each crate README
for the MPI status of bond/family construction.

## License

MIT OR Apache-2.0
