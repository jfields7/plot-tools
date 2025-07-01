# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2025, Jacob Fields <jmf6719@psu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from .load_ath_bin import *
from scipy.interpolate import griddata

class Image:
  """
    This class represents an image constructed from a BinaryData file. Effectively what it
    does is it creates an single interpolator from all the MeshBlock data in a specified
    region in a feeble attempt to prevent the renderer from doing weird things because a
    plot is broken up into hundreds of small pieces.
    TODO: Right now this assumes 2D data.
  """
  def __init__(self, filename, extent, var):
    self.data = BinaryData(filename)
    self.extent = extent
    self.var = var

    self.ConstructData()

  def ConstructData(self):
    xpoints = np.empty(0)
    ypoints = np.empty(0)
    data    = np.empty(0)

    with open(self.data.filename, 'rb') as f:
      for block in self.data.blocks:
        xs, ys = block.get_coord_blocks()
        block_extent = block.get_extent()
        # Exclude blocks that will not be in the final image.
        if self.extent[1] < block_extent[0] or self.extent[0] > block_extent[1] \
           or self.extent[3] < block_extent[2] or self.extent[2] > block_extent[3]:
          continue

        xpoints = np.append(xpoints, xs.flatten(order='C'))
        ypoints = np.append(ypoints, ys.flatten(order='C'))

        block_data = block.get_var(f, self.var, self.data.block_cell_format).transpose()
        data = np.append(data, block_data.flatten(order='C'))
    
    self.xpoints = xpoints
    self.ypoints = ypoints
    self.idata = data

  def transform_coordinates(self, xtransform, ytransform):
    self.extent[0:2] = xtransform([self.extent[0], self.extent[1]])
    self.extent[2:4] = ytransform([self.extent[2], self.extent[3]])

    self.xpoints = xtransform(self.xpoints)
    self.ypoints = ytransform(self.ypoints)

  def make_image_data(self, xres, yres, method='nearest', return_grid=False,
                      transform=None):
    # Make the interpolation grid
    xcoord = np.linspace(self.extent[0], self.extent[1], xres)
    ycoord = np.linspace(self.extent[2], self.extent[3], yres)
    X, Y = np.meshgrid(xcoord, ycoord, indexing='ij')
    xpoints = self.xpoints
    ypoints = self.ypoints
    if not transform == None:
      ypoints, xpoints = transform(ypoints, xpoints)
    if return_grid:
      return X, Y, griddata((xpoints, ypoints), self.idata, (X, Y),
                            method=method).transpose()
    return griddata((xpoints, ypoints), self.idata, (X, Y),
                    method=method).transpose()
