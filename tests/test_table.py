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
import pytest

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.utils.eos import EquationOfState

def test_load_logtable():
  try:
    table = EquationOfState('tables/SFHo_reduced.h5')
  except Exception as e:
    pytest.fail(f'Unexpected {err=}, {type(err)=}')

def test_load_NQTtable():
  try:
    table = EquationOfState('tables/SFHo_reduced_NQT.h5')
  except Exception as e:
    pytest.fail(f'Unexpected {err=}, {type(err)=}')
