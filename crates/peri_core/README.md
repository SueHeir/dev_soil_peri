# peri_core

The batteries-included umbrella crate for [dev_soil_peri](../../README.md) — depend on
this.

```rust
use peri_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // input, comm, domain, neighbour, run, output
       .add_plugins(BondBasedDefaultPlugins); // point cloud, Verlet, bond-based PD force
    app.start();
}
```

Re-exports every dev_soil_peri plugin crate (`peri_atom`, `peri_bond`) and the SOIL
substrate crates it builds on, and provides two plugin groups:

| Group | Purpose |
|---|---|
| `CorePlugins` | infrastructure: input, communication, domain decomposition, neighbour lists, groups, run loop, output |
| `BondBasedDefaultPlugins` | bond-based PD: `PeriAtom` + materials, lattice point cloud, velocity-Verlet integration, nonlocal bond force with fracture |

The `prelude` re-exports everything a typical simulation needs.
