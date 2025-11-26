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
from matplotlib.colors import LogNorm

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.load_sph_vtk import SphericalData

data = SphericalData('data/bns_g2.r=25.00.mhd_w_bcc.00438.vtk')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surf, sm = data.plot_data('dens', ax, cmap='inferno', norm=LogNorm(vmin=3e-15, vmax=1.28e-3))

ax.set_title(f't = {data.time}')

fig.colorbar(sm, ax=ax)

plt.show()
