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

import matplotlib.pyplot as plt

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.tracker import Tracker

co_0 = Tracker('data/bns.co_0.txt')
co_1 = Tracker('data/bns.co_1.txt')

x0, y0 = co_0.smooth_walk()
x1, y1 = co_1.smooth_walk()

plt.plot(co_0.x, co_0.y, 'm.', label='BNS 1, unsmoothed')
plt.plot(co_1.x, co_1.y, 'c.', label='BNS 2, unsmoothed')
plt.plot(x0, y0, 'r-', label='BNS 1, smoothed')
plt.plot(x1, y1, 'b-', label='BNS 2, smoothed')

plt.xlabel(r'$x~[\mathrm{M}_\odot]$')
plt.ylabel(r'$y~[\mathrm{M}_\odot]$')
plt.gca().set_aspect('equal')
plt.legend()

plt.show()
