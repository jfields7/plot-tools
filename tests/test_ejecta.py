import numpy as np

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.utils.eos import EquationOfState
from athplot.utils.ejecta import Ejecta

stem = 'data/bns.r=400.00.'
end  = '.01640.vtk'

eos = EquationOfState('tables/SFHo_reduced_NQT.h5', use_NQT=True)
ejecta = Ejecta(stem + 'mhd_w_bcc' + end, stem + 'adm' + end, eos, fname_gauge=(stem + 'z4c' + end))

Mej_rate, Mej_rate_Ber, Mej_rate_geo = ejecta.calc_Mej_rate()

print(Mej_rate)
print(Mej_rate_Ber)
print(Mej_rate_geo)
