# CHANGELOG

# 0.3

## New Features and Enhancements
* Spherical VTK data now reads comment metadata
* Trackers now smooth using the midpoint rather than the first point, which greatly
  reduces the sawtooth artifacts in walking trackers.

## Bug Fixes
* Fixed a typo in `get_coord_blocks` for sliced binary data.
* `plot_slice` for `BinaryData` and `GroupData` no longer generates blocking artifacts
  for logarithmically-scaled colorbars when the dynamical range spans 10+ orders of
  magnitude.
