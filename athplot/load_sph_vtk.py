# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2025, Jacob Fields <fields@ias.edu>
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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from vtk.util import numpy_support
from vtk import vtkStructuredGridReader

class SphericalData:
  """
    This class represents a VTK spherical data dump from AthenaK.

    A dump may contain several concentric surfaces sharing one angular grid
    (an <output> block with 'radii = ...' or 'r_min/r_max/nradii'), or a single
    surface ('radius = ...'), which is just the nradii == 1 case. Files written
    before multiple radii were supported are read as well.

    Note that this utility relies on the VTK python package.
  """
  def __init__(self, filename, encoding='latin_1'):
    # Extract metadata stored as a comment in the second line
    string = ""
    with open(filename, 'r', encoding=encoding) as f:
        l = 0
        for line in f:
            if l == 1:
                string = line.strip()
                break
            else:
                l += 1

    # Defaults for dumps written before multiple radii were supported
    self.nradii = 1
    self.rmin = None
    self.rmax = None

    substrings = string.split()
    for s in substrings:
        if s[:5] == 'time=':
            self.time = float(s[5:])
        if s[:6] == 'cycle=':
            self.cycle = int(s[6:])
        if s[:7] == 'nradii=':
            self.nradii = int(s[7:])
        if s[:5] == 'rmin=':
            self.rmin = float(s[5:])
        if s[:5] == 'rmax=':
            self.rmax = float(s[5:])
        if s[:4] == 'rad=':
            self.rmin = float(s[4:])
            self.rmax = self.rmin
        if s[:3] == 'xc=':
            self.xc = float(s[3:])
        if s[:3] == 'yc=':
            self.yc = float(s[3:])
        if s[:3] == 'zc=':
            self.zc = float(s[3:])

    # Open the file with the VTK reader
    reader = vtkStructuredGridReader()

    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()

    reader.update()

    # Read in the coordinate data
    cells_vtk = reader.GetOutput().GetPoints().GetData()
    cells = numpy_support.vtk_to_numpy(cells_vtk)

    # Reorganize the data into separate r, theta, and phi quantities
    dims = [0, 0, 0]
    reader.GetOutput().GetDimensions(dims)
    self.nradii = dims[0]
    self.ntheta = dims[1]
    self.nphi = dims[2]

    # Shape of a single shell, and of a full field including the radial axis
    self.shape = (dims[2], dims[1])
    self.shape3d = (dims[0], dims[2], dims[1])

    # Every point stores (radius, theta, phi) and the radius varies fastest
    self.radii = cells[:self.nradii,0]
    self.theta = np.reshape(cells[::self.nradii,1], self.shape)
    self.phi = np.reshape(cells[::self.nradii,2], self.shape)

    if self.rmin is None:
      self.rmin = self.radii[0]
    if self.rmax is None:
      self.rmax = self.radii[-1]

    # Read in the field data
    nfields = reader.GetOutput().GetPointData().GetNumberOfArrays()

    self.field_names = []
    self.fields3d = {}
    for i in range(nfields):
      field = reader.GetOutput().GetPointData().GetArrayName(i)
      self.field_names.append(field)
      field_data = reader.GetOutput().GetPointData().GetArray(i)
      self.fields3d[field] = np.moveaxis(np.reshape(numpy_support.vtk_to_numpy(field_data),
                                                    self.shape + (self.nradii,)), -1, 0)

    # Distinguish between nradii == 1 and nradii > 1 for backward-compatibility
    if self.nradii == 1:
      self.r = self.radii[0]
      self.fields = {k: v[0] for k, v in self.fields3d.items()}
    else:
      self.r = self.radii
      self.fields = self.fields3d

  # Return var on a single shell, with shape (nphi, ntheta).
  def shell(self, var, iradius=0):
    return self.fields3d[var][iradius]

  # Reduce the var over the angles of each shell, giving radial profile (nradii,)
  def radial_profile(self, var, reduction=np.mean):
    return reduction(self.fields3d[var].reshape(self.nradii, -1).astype(np.float64),
                     axis=1)

  # Integrate var over each surface, giving array of (nradii,)
  def integrate(self, var):
    if 'weights' not in self.fields3d:
      raise RuntimeError("This dump does not contain a 'weights' array.")
    w = self.fields3d['weights'].reshape(self.nradii, -1)
    integrand = (w*self.fields3d[var].reshape(self.nradii, -1)).astype(np.float64)
    return np.sum(integrand, axis=1)

  def plot_data(self, var, ax, cmap='viridis', norm=None, vmin=None, vmax=None,
                iradius=0):
    # Note that phi is periodic, so for plotting purposes we need to append an extra
    # column of data to get a closed surface
    r = self.radii[iradius]
    theta = np.append(self.theta, [self.theta[0]], axis=0)
    phi = np.append(self.phi, [2*np.pi*np.ones(self.shape[1])], axis=0)

    data = np.copy(self.fields3d[var][iradius])
    data = np.append(data, [data[0]], axis=0)

    X = r*np.sin(theta)*np.cos(phi)
    Y = r*np.sin(theta)*np.sin(phi)
    Z = r*np.cos(theta)

    if norm != None:
      if vmin != None:
        raise RuntimeError("Cannot specify vmin and norm at the same time!")
      if vmax != None:
        raise RuntimeError("Cannot specify vmax and norm at the same time!")
      vmin = norm.vmin
      vmax = norm.vmax

    if vmin==None:
      vmin = data.min()
    if vmax==None:
      vmax = data.max()

    if norm == None:
      norm = Normalize(vmin=vmin, vmax=vmax)

    plot_data = norm(data)

    colors = plt.get_cmap(cmap)(plot_data)

    surf = ax.plot_surface(X, Y, Z, shade=False, facecolors=colors)

    sm = ScalarMappable(norm=norm, cmap=cmap)

    return surf, sm
