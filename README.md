# plot-tools
A collection of plotting tools and utilities for AthenaK. Here are the files currently
included:
* `load_ath_bin.py`: Defines the `BinaryData` and `MeshBlock` classes. Contains
  functionality for plotting slices of data, derived variables, etc.
* `image.py`: Defines the `Image` class, which can interpolate a 2D slice of `BinaryData`
  onto a Cartesian grid for simpler manipulation. Generally slower than `BinaryData`.
* `tracker.py`: Defines the `Tracker` class, which can read and manipulate the data output
  by AthenaK's compact object trackers.

Here are some limitations of the tools:
* `BinaryData` streams data off the disk a mesh block at a time. While this is more
  efficient for plotting very large datasets, it means that data ranges must be specified
  beforehand. Additionally, limitations in `matplotlib` mean that plot ranges spanning
  many orders of magnitude are typically truncated to the data range of the each block,
  which can result in artifacts.
* Many features, such as derived variables and `Image`, are only supported for 2D data.
* The `smooth_walk` function of Tracker currently assumes the trajectory of interest is
  only in the x-y plane.

A number of examples are provided in the `examples` directory.

Finally, a disclaimer: these tools are not polished and may be quite clumsy to use. Use
at your own risk.
