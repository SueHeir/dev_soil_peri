# dev_soil_peri

<!-- disclaimer-banner -->
> This code was fully written via **Claude 4.6,4.8, Fable 5, and Codex**, and stands as a proof of concept for a **bevy-like** ecosystem for physics simulation research, with the goal of testing if one scheduler/framework (**GRASS**) works for most scientific codes. **SOIL** and **FIELD** are particle- and mesh-based substrates for physics such as **DIRT** (DEM) or **dev_field_efvm**. Note that all other physics based repos I have start with **dev_**, as I do **NOT** know these methods. Please read, evaluate, use with a grain of salt, I have not personally read or reviewed everything here.
<!-- /disclaimer-banner -->


Peridynamics (PD) physics in the **GRASS → SOIL → {DIRT, dev_soil_peri}** stack.

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists
       ├─ DIRT            DEM physics: contact, bonds, walls, clumps
       └─ dev_soil_peri   peridynamics: nonlocal bond forces, fracture, damage   ← you are here
```

dev_soil_peri is the peridynamics sibling of [DIRT](https://github.com/SueHeir/dirt). Where
DIRT resolves short-range *contact* between discrete grains, dev_soil_peri resolves the
*nonlocal* interaction of continuum material points: every point interacts with
every neighbour inside a finite **horizon** δ through a bond, and material
failure emerges naturally as those bonds break. It rides the
[SOIL](https://github.com/SueHeir/soil) substrate for all the method-agnostic
machinery — atom data, domain decomposition, halo exchange, neighbour lists,
time integration — and adds only the peridynamic constitutive models.

The fit is close to exact: a PD point cloud *is* a SOIL atom set, the PD horizon
*is* a SOIL neighbour cutoff, and the PD family (the fixed set of bonds formed in
the reference configuration) *is* a SOIL [`BondStore`]. dev_soil_peri therefore reuses
SOIL's persistent, tag-addressed, MPI-aware bond storage rather than inventing
its own.

## Composition boundary

Peridynamic resources and bond semantics belong here; SOIL owns their generic
particle transport, and GRASS holds the resource instances and schedules the
plugins. GRASS does not define a material point or bond. If this solver is
coupled later, the coupling belongs in a separate parent plugin which drives
the child to exported scheduler seams and exchanges only the public interface
state. No child-level cross-App `MultiRes` access is required or supported.

This repository currently validates its single-solver, single-process slice. It
does not yet claim a topology-driven multi-rank coupling example.

## Layout

`peri_core` is the batteries-included umbrella crate (depend on this):

```rust
use peri_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
       .add_plugins(BondBasedDefaultPlugins);
    app.start();
}
```

| crate | role |
|---|---|
| [`peri_core`](crates/peri_core/README.md) | umbrella: `CorePlugins`, `BondBasedDefaultPlugins`, prelude |
| [`peri_atom`](crates/peri_atom/README.md) | per-point PD data (`PeriAtom`: reference position, volume, damage), PD materials, lattice point-cloud generation |
| [`peri_bond`](crates/peri_bond/README.md) | bond-based peridynamics: constant-micromodulus pairwise force, critical-stretch bond breakage, damage field |

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
> the ν = 1/4 restriction, is the planned next tier (`peri_state`) on this same
> substrate.

## Building

Clone dev_soil_peri and build — its [soil](https://github.com/SueHeir/soil) and
[grass](https://github.com/SueHeir/grass) dependencies are pulled from GitHub
automatically.

```bash
git clone https://github.com/SueHeir/dev_soil_peri
cd dev_soil_peri
cargo run --release --example bb_bar_tension --no-default-features -- examples/bb_bar_tension/config.toml
```

The `--no-default-features` flag runs single-process (no MPI), which is the
supported configuration for the current bond-based slice — see each crate README
for the MPI status of bond/family construction.

## Examples

| example | what it shows |
|---|---|
| [`bb_bar_tension`](examples/bb_bar_tension) | elastic bond-based bar in uniaxial tension (nonlocal stress wave, no failure) |
| [`bb_plate_fracture`](examples/bb_plate_fracture) | **mode-I fracture of a pre-notched plate** — a crack nucleates at the notch tip and propagates straight along the symmetry plane by critical-stretch bond breakage (nothing prescribes the path). Validates the Griffith relation `s₀ = √(5G₀/9Kδ)` and an emergent propagation threshold. |

```bash
# mode-I pre-notched plate fracture (single-process)
cargo run --release --example bb_plate_fracture --no-default-features -- examples/bb_plate_fracture/config.toml
# end-to-end validation (build + run + checks + validation plot)
python3 examples/bb_plate_fracture/verify.py
```

## License

MIT OR Apache-2.0
