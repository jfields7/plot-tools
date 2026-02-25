# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2026, Jacob Fields <fields@ias.edu>
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
import h5py as h5

from scipy.interpolate import RegularGridInterpolator, interp1d

class EquationOfState:
  """
    This class represents a PyCompOSE equation of state.
    Inputs:
      filename: the path to the PyCompOSE HDF5 table.
      use_NQT: whether or not the table is spaced in second-order NQT space (see
               https://arxiv.org/abs/2501.05410). Default: False
  """
  def __init__(self, filename, use_NQT=False):
    self.log2_ = np.log2
    self.exp2_ = np.exp2
    if use_NQT:
      from compose.NQTs.NQTLib import NQT_exp2_O2 as NQT_exp
      from compose.NQTs.NQTLib import NQT_log2_O2 as NQT_log

      self.log2_ = NQT_log
      self.exp2_ = NQT_exp

    table = h5.File(filename,'r')

    # Scalars
    self.mn = table['mn'][()]

    # Independent variables
    self.lnb = self.log2_(table['nb'][:])
    self.yq  = table['yq'][:]
    self.lt  = self.log2_(table['t'][:])

    # Input limits
    self.nb_min = table['nb'][0]
    self.nb_max = table['nb'][-1]
    self.yq_min = table['yq'][0]
    self.yq_max = table['yq'][-1]
    self.t_min  = table['t'][0]
    self.t_max  = table['t'][-1]

    # Important fields
    press = np.array(table['Q1'])*(np.array(table['nb'])[:,np.newaxis,np.newaxis])
    e = (np.array(table['Q7']) + 1.)*(np.array(table['nb'])[:,np.newaxis,np.newaxis])*self.mn
    cs2 = table['cs2']

    self.press_raw = press
    self.e_raw = e
    self.cs2_raw = cs2

    # Derived scalars
    self.hmin = np.min((e + press)/(np.array(table['nb'])[:,np.newaxis,np.newaxis]))

    # Construct interpolators for fields
    self.lpress = RegularGridInterpolator((self.lnb, self.yq, self.lt), self.log2_(press),
                   method='linear', bounds_error=False, fill_value=None)
    self.le     = RegularGridInterpolator((self.lnb, self.yq, self.lt), self.log2_(e),
                   method='linear', bounds_error=False, fill_value=None)
    self.csq   = RegularGridInterpolator((self.lnb, self.yq, self.lt), cs2,
                   method='linear', bounds_error=False, fill_value=None)
  
  def calc_temp_from_press(self, nb, yq, press):
    pass

  def calc_temp_from_e(self, nb, yq, e):
    pass

  """
    Compute the specific enthalpy from the number density, charge fraction, and
    temperature.

    Inputs:
      nb: an array of number density points
      yq: an array of charge fraction points
      t:  an array of temperature points

    All three inputs are expected to have the same shape.

    Returns:
      an array of the same shape as the inputs.
  """
  def calc_enthalpy(self, nb, yq, t):
    pts = np.column_stack((
      self.log2_(np.clip(nb.flatten(), self.nb_min, self.nb_max)),
      np.clip(yq.flatten(), self.yq_min, self.yq_max),
      self.log2_(np.clip(t.flatten(), self.t_min, self.t_max))))

    p = self.exp2_(self.lpress(pts))
    e = self.exp2_(self.le(pts))
    return ((e + p)/(np.clip(nb.flatten(), nb_min, nb_max)*self.mn)).reshape(nb.shape)

  """
    Compute the pressure from the number density, charge fraction, and
    temperature.

    Inputs:
      nb: a Numpy array of number density points
      yq: a Numpy array of charge fraction points
      t:  a Numpy array of temperature points

    All three inputs are expected to have the same shape.

    Returns: a Numpy array of the same shape as the inputs.
  """
  def calc_pressure(self, nb, yq, t):
    pts = np.column_stack((
      self.log2_(np.clip(nb.flatten(), self.nb_min, self.nb_max)),
      np.clip(yq.flatten(), self.yq_min, self.yq_max),
      self.log2_(np.clip(t.flatten(), self.t_min, self.t_max))))

    return self.exp2_(self.lpress(pts)).reshape(nb.shape)

  """
    Compute the total energy density from the number density, charge fraction, and
    temperature.

    Inputs:
      nb: a Numpy array of number density points
      yq: a Numpy array of charge fraction points
      t:  a Numpy array of temperature points

    All three inputs are expected to have the same shape.

    Returns: a Numpy array of the same shape as the inputs.
  """
  def calc_energy_dens(self, nb, yq, t):
    pts = np.column_stack((
      self.log2_(np.clip(nb.flatten(), self.nb_min, self.nb_max)),
      np.clip(yq.flatten(), self.yq_min, self.yq_max),
      self.log2_(np.clip(t.flatten(), self.t_min, self.t_max))))

    return self.exp2_(self.le(pts)).reshape(nb.shape)
