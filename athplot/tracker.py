# AthPlot: manipulate and plot AthenaK binary data                                                     # Copyright (C) 2025, Jacob Fields <jmf6719@psu.ed>                                                    #                                                                                                      # This program is free software: you can redistribute it and/or modify                                 # it under the terms of the GNU General Public License as published by                                 # the Free Software Foundation, either version 3 of the License, or                                    # (at your option) any later version.                                                                  #                                                                                                      # This program is distributed in the hope that it will be useful,                                      # but WITHOUT ANY WARRANTY; without even the implied warranty of                                       # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                        # GNU General Public License for more details.                                                         #                                                                                                      # You should have received a copy of the GNU General Public License                                    # along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy.interpolate import make_smoothing_spline

class Tracker:
  def __init__(self, fname):
    # Load the data, carefully ensuring that duplicate data is removed
    data = np.genfromtxt(fname)
    self.time, idcs = np.unique(data[:,1], return_index=True)
    self.it = data[idcs,0]
    self.x = data[idcs,2]
    self.y = data[idcs,3]
    self.z = data[idcs,4]
    self.vx = data[idcs,5]
    self.vy = data[idcs,6]
    self.vz = data[idcs,7]

  # Construct interpolating splines to the trajectories created by
  # 'walk' mode trackers.
  def smooth_walk(self):
    # First, we need to clean up the data for the x and y directions, which
    # we do separately since they often change at different times
    x_idcs = [0]
    for i in range(1,len(self.time)):
      if not self.x[x_idcs[-1]] == self.x[i]:
        x_idcs.append(i)

    y_idcs = [0]
    for i in range(1,len(self.time)):
      if not self.y[y_idcs[-1]] == self.y[i]:
        y_idcs.append(i)

    # Now construct interpolating splines for each trajectory
    cs_x = make_smoothing_spline(self.time[x_idcs], self.x[x_idcs])
    cs_y = make_smoothing_spline(self.time[y_idcs], self.y[y_idcs])

    return cs_x(self.time), cs_y(self.time)
