# peri_bond

Bond-based peridynamics for [dev_soil_peri](../../README.md).

Each material point bonds, once in the reference configuration, to every point
inside the horizon δ. Bonds are stored in SOIL's `BondStore` (partner addressed
by global tag). Every step each point sums the pairwise bond force on itself:

```text
F_i = V_i · Σ_j  c · s · (y_j − y_i)/|y_j − y_i| · V_j
```

with stretch `s = (|y| − r₀)/r₀` and constant micromodulus `c`. A bond breaks
irreversibly once `s` exceeds the critical stretch `s₀`; the fraction of a
point's broken bonds is its **damage**, and a crack is the locus where damage → 1.

Each point applies the force only to itself (a bond is evaluated by both
endpoints), so the force loop needs no Newton's-third-law reverse communication.

**Status:** family construction walks SOIL's binned neighbour/halo data in the
reference configuration. Families may include cross-rank partners, addressed by
global tag, while bond topology stays in SOIL's migrating/restartable
`BondStore`.
