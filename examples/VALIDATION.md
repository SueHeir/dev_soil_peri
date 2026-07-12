# Validation Examples

- [bb_plate_fracture](bb_plate_fracture/) checks the bond-based fracture gate
  against the Silling-Askari critical-stretch relation and a digitized
  Ha-Bobaru crack trajectory. Its `verify_decomposition.py` independently
  reruns that serial gate before comparing the same fracture case at one and
  two MPI ranks; scope and limits are documented beside the example.
  Its [`decomposition_invariance.svg`](bb_plate_fracture/plots/decomposition_invariance.svg)
  is generated only after both gates pass and shows the observed MPI deltas
  against the stated f64 reduction bounds.
- [peri_family_consistency](peri_family_consistency/) checks one-rank vs
  two-rank SOIL neighbor/halo-backed family construction, damage, and crack
  span.
