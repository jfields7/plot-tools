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

# A collection of functions for integrating over spherical surfaces. Derived from scripts
# from Ruocheng Zhai.

import numpy as np

from ..load_sph_vtk import SphericalData

def calc_ur(velx, vely, velz, theta, phi):
  ur = (velx * np.sin(theta) * np.cos(phi) +
        vely * np.sin(theta) * np.sin(phi) +
        velz * np.cos(theta))

  return ur

def sqrt_det_metric_adm(gxx, gxy, gxz, gyy, gyz, gzz, alpha):
  gamma_det = (gxx * gyy * gzz -
               2 * gxy * gxz * gyz -
               gxx * gyz**2 -
               gyy * gxz**2 -
               gzz * gxy**2)
  return np.sqrt(gamma_det) * alpha

def integrate_over_sphere(integrand, theta_1d, phi_1d):
  def get_widths(c):
    # Calculate cell boundaries from centers
    e = np.zeros(len(c) + 1)
    e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0] = c[0] - (e[1] - e[0])
    e[-1] = c[-1] + (c[-1] - e[-2])
    return np.abs(np.diff(e)) # Abs handles decreasing coordinate sequences

  dt = get_widths(theta_1d) # Shape (Ntheta,)
  dp = get_widths(phi_1d)   # Shape (Nphi,)

  # Create 2D weight map for poles
  # We use broadcasting to check theta values across the whole grid
  weights = np.ones_like(integrand)
  is_pole = np.isclose(theta_1d, 0, atol=1e-7) | np.isclose(theta_1d, np.pi, atol=1e-7)

  # weights is (Nphi, Ntheta), is_pole is (Ntheta,)
  # Divide the columns corresponding to poles by Nphi
  weights[:, is_pole] /= len(phi_1d)

  # Riemann sum using broadcasting: (Nphi, Ntheta) * (1, Ntheta) * (Nphi, 1)
  return np.sum(integrand * weights * dt[None, :] * dp[:, None])

def integrate_over_phi(integrand, theta_1d, phi_1d):
  def get_widths(c):
    # Calculate cell boundaries from centers
    e = np.zeros(len(c) + 1)
    e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0] = c[0] - (e[1] - e[0])
    e[-1] = c[-1] + (c[-1] - e[-2])
    return np.abs(np.diff(e)) # Abs handles decreasing coordinate sequences

  dp = get_widths(phi_1d) # Shape (Nphi,)

  # Create 2D weight map for poles
  # integrand shape is typically (Nphi, Ntheta)
  weights = np.ones_like(integrand)
  is_pole = np.isclose(theta_1d, 0, atol=1e-7) | np.isclose(theta_1d, np.pi, atol=1e-7)

  # Divide the columns corresponding to poles by Nphi
  # (Since we sum over phi, we only average the contribution at the singular poles)
  weights[:, is_pole] /= len(phi_1d)

  # Perform the sum over the phi axis (axis 0)
  # result shape: (Ntheta,)
  return np.sum(integrand * weights * dp[:, None], axis=0)
