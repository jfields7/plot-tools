import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(SCRIPTDIR, os.pardir))
from athplot.load_ath_bin import BinaryData

data = BinaryData("data/bns.mhd_w_bcc.00100.bin")

pcm = data.plot_slice("dens", cmap='inferno', norm=LogNorm(vmin=1.28e-13, vmax=1.28e-3),
                      interpolation='nearest')

plt.colorbar(pcm, label=r'$\rho$')

plt.xlim([-80, 80])
plt.ylim([-80, 80])
plt.title(f'Time: {data.time:.2f} Msun')
plt.show()
