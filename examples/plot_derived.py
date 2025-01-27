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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.load_ath_bin import BinaryData

def calc_lorentz(ux, uy, uz):
  return np.sqrt(1 + ux*ux + uy*uy + uz*uz)

data = BinaryData('data/blast_table.mhd_w_bcc.00200.bin')

data.register_derived_variable("derived:lorentz", calc_lorentz, "velx", "vely", "velz")

pcm = data.plot_slice("derived:lorentz", cmap='plasma',
                      norm=Normalize(vmin=1.0, vmax=1.6),
                      interpolation='nearest')

plt.colorbar(pcm, label=r'$W$', extend='max')

plt.xlim([-6, 6])
plt.ylim([-6, 6])
plt.show()
