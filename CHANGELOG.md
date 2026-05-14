# CHANGELOG

# 0.3

## New Features and Enhancements
* Spherical VTK data now reads comment metadata
* Trackers now smooth using the midpoint rather than the first point, which greatly
  reduces the sawtooth artifacts in walking trackers.
* New EOS tooling for reading
  [PyCompOSE](https://github.com/computationalrelativity/PyCompOSE) tables
* Automated tests for EOS
* New `Ejecta` class for computing ejecta properties from spherical VTK data
* Added `make_var_data` function to `Image` to simplify reading multiple variables
  from the same `Image`.
* New `plot_line` function in `Image` class to interpolate data along a specified
  line.
* `BinaryData` constructor now also reads in cycle count.

## Bug Fixes
* Fixed a typo in `get_coord_blocks` for sliced binary data.
* `plot_slice` for `BinaryData` and `GroupData` no longer generates blocking artifacts
  for logarithmically-scaled colorbars when the dynamical range spans 10+ orders of
  magnitude.
* Fixed issue with axis ordering in `Image`
