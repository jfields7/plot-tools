# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2025, Jacob Fields <jmf6719@psu.ed>
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

class GroupData:
  """
    This class contains a collection of BinaryData objects which all correspond to the
    same dataset at the same time. This is designed to make it easier to manipulate
    datasets with variables spanning multiple output files. For example, DynGRMHD datasets
    will store the metric in ADM output and the magnetic field in a primitive or conserved
    output, so computing b^2/rho will require both the primitive and ADM outputs.
  """
  def __init__(self, filenames):
    # Collect all the datasets
    self.datasets = []
    for f in filenames:
      self.datasets.append(BinaryData(f))

    # Validate that the datasets are consistent with one another
    nblocks = len(self.datasets[0].blocks)
    size = self.datasets[0].blocks[0].size
    is_2d = self.datasets[0].is_2d()
    for data in self.datasets:
      ds_nblocks = len(data.blocks)
      if ds_nblocks != nblocks:
        raise RuntimeError('Not all datasets have the same number of meshblocks.')

      ds_size = data.blocks[0].size
      if ds_size != size:
        raise RuntimeError('Not all datasets have the same size meshblocks.')

      ds_is_2d = data.is_2d()
      if ds_is_2d != is_2d:
        raise RuntimeError('Not all datasets have the same dimensionality.')

    self.derived_vars = {}
    
  def register_derived_variable(self, name, f, *args):
    if len(name) < 8:
      raise RuntimeError("Derived variable name must begin with 'derived:'!")
    elif not name[:8] == 'derived:':
      raise RuntimeError("Derived variable name must begin with 'derived:'!")
    self.derived_vars[name] = [f, *args]

  def plot_slice(self, var, ax=None, cmap='viridis', norm=None, vmin=None, vmax=None,
                 interpolation='none', origin='lower', slice_loc=None, rescale=1.):
    pcm = None
    if ax == None:
      ax = plt.gca()

    if len(var) >= 8:
      if var[:8] == 'derived:':
        if self.datasets[0].is_2d():
          # Determine which variables are in which datasets
          fn = self.derived_vars[var][0]
          fnvars = self.derived_vars[var][1:]
          dsets = []
          for var in fnvars:
            for i, data in enumerate(self.datasets):
              if var in data.variable_names:
                # Variable is in dataset, stop searching for this one
                dsets.append(i)
                break
              else:
                # Keep searching this dataset
                continue
          
          # Now open all the datasets at once
          files = []
          for data in self.datasets:
            files.append(open(data.filename, 'rb'))

          # Iterate over all the blocks
          for i in range(0,len(self.datasets[0].blocks)):
            fnargs = []
            # Extract each variable from its corresponding dataset and block
            for v,d in zip(fnvars, dsets):
              fnargs.append(self.datasets[d].blocks[i].get_var(files[d], v,
                                        self.datasets[d].block_cell_format)*rescale)
            
            # Calculate and plot the derived variable
            img_data = fn(*tuple(fnargs))
            extent = self.datasets[0].blocks[i].get_extent()
            pcm = ax.imshow(img_data, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
                       interpolation=interpolation, origin=origin,
                       extent=extent, interpolation_stage='rgba')

          for f in files:
            f.close()
        else:
          raise RuntimeError('Derived variables not currently supported for 3D slicing.')
        return pcm
    
    # If we aren't plotting a derived variable, we can simply resort to normal plotting
    for data in self.datasets:
      if var in data.variable_names:
        pcm = data.plot_slice(var, ax, cmap, norm, vmin, vmax, interpolation, origin,
                              slice_loc, rescale)
        break
      else:
        continue
    
    return pcm
