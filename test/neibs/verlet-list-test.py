from mdcraft.tools import Threads
threads = Threads(0)
from mdcraft.data import Atoms
from mdcraft.neibs import VerletList

import numpy as np
import math
import matplotlib.pyplot as plt
import time as Time

n = 1000 * 8;
volume = 1. * 1.
atoms = Atoms(5)
neibs     = Atoms(n)
d = math.pow(volume / n, 1. / 2.)

np.random.seed(n)
print("Creating random ", n, " points in storage neibs...  ")

neibs.data["uid"] = 0
atoms.data["r"]["x"] = [-0.4, -0.2, 0.0, 0.2, 0.4]
atoms.data["r"]["x"] -= 0.08
atoms.data["r"]["y"] = 0.0
atoms.data["r"]["z"] = 0.0

neibs.data["uid"] = np.arange(1, n + 1)
neibs.data["r"]["x"] = np.random.rand(n) - 0.5
neibs.data["r"]["y"] = np.random.rand(n) - 0.5
atoms.data["r"]["z"] = 0.0
r = 4.0*d
neibs.data["rns"] = r

np.random.seed(int(Time.time()))
print("Constructing neibslist...  ")

nlist = VerletList(
	atoms      = atoms,
	neighbors  = neibs,
	periodic_x = True, 
	periodic_y = True,
	threads    = threads
)

start = Time.time()

nlist.update()

elapsed = Time.time() - start
print("elapsed time: ", elapsed)


plt.scatter(
	neibs.data["r"]["x"], 
	neibs.data["r"]["y"], 
	color="blue", 
	s = 40
)
ax = plt.subplot()
for i in range(atoms.data.size):
	ineibs = nlist[i]
	print(ineibs, ineibs.size)
	plt.scatter(
		neibs.data["r"]["x"][ineibs], 
		neibs.data["r"]["y"][ineibs], 
		color="red", 
		s= 40
	)
	x = atoms.data["r"]["x"][i]
	y = atoms.data["r"]["y"][i]
	plt.scatter(x, y, color="black", s = 80)
	c0 = plt.Circle((x, y), r, fill = False, color = 'b')
	ax.add_artist(c0)

plt.axis("equal")
plt.grid()
plt.tight_layout()
plt.show()
