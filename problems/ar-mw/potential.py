import numpy as np
import matplotlib.pyplot as plt

from mdcraft.tools import Threads
from mdcraft.solver.potential import LJs

Rcutoff = 0.8125

potential = LJs(
	aVr     = 1.03120074782442750e-3, # kJ/mol
	rVr     = 0.33841043857528683, 
	Rcutoff = 0.8125 # nm
)

r = np.linspace(0.3, Rcutoff, 100)
r2 = r*r

plt.plot(r, potential.value(r*r))
plt.grid(True)
plt.show()