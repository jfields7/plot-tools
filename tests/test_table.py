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

compose = pytest.importorskip('compose.utils')

def test_load_logtable():
  try:
    table = EquationOfState('tables/SFHo_reduced.h5')
  except Exception as e:
    pytest.fail(f'Unexpected {err=}, {type(err)=}')


def test_load_NQTtable():
  assert compose is not None
  try:
    table = EquationOfState('tables/SFHo_reduced_NQT.h5',use_NQT=True)
  except Exception as e:
    pytest.fail(f'Unexpected {err=}, {type(err)=}')

def check_pressure(use_NQT):
  if use_NQT:
    table = EquationOfState('tables/SFHo_reduced_NQT.h5',use_NQT=True)
  else:
    table = EquationOfState('tables/SFHo_reduced.h5')
  lnb, yq, lt = np.meshgrid(table.lnb, table.yq, table.lt, indexing='ij')

  p_calc = table.calc_pressure(table.exp2_(lnb), yq, table.exp2_(lt))

  err = np.abs(p_calc - table.press_raw)
  return np.max(err) == pytest.approx(0.0, rel=1e-15)

def check_energy(use_NQT):
  if use_NQT:
    table = EquationOfState('tables/SFHo_reduced_NQT.h5',use_NQT=True)
  else:
    table = EquationOfState('tables/SFHo_reduced.h5')
  lnb, yq, lt = np.meshgrid(table.lnb, table.yq, table.lt, indexing='ij')

  e_calc = table.calc_energy_dens(table.exp2_(lnb), yq, table.exp2_(lt))

  err = np.abs(e_calc - table.e_raw)
  return np.max(err) == pytest.approx(0.0, rel=1e-15)

def test_logpressure():
  assert check_pressure(False) == True

def test_NQTpressure():
  assert compose is not None
  assert check_pressure(True) == True

def test_logenergy():
  assert check_energy(False) == True

def test_NQTenergy():
  assert compose is not None
  assert check_energy(True) == True
