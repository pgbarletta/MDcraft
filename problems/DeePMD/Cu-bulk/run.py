import numpy as np
import math, sys, os, time as Time
import h5py as hdf
from mdcraft.tools import Threads
threads = Threads(4)

from mdcraft.data import Atoms
from mdcraft.neibs import VerletList
from mdcraft.lattice import Domain

from mdcraft.solver            import ManyBodySolver as Solver
from mdcraft.solver.potential  import DeepModelingPotential as Potential

from mdcraft.solver.thermostat import Langevin
from mdcraft.solver.stepper    import Verlet as Stepper

kBuf	 = 1.25           # buffer size
Ni_mass  = 58.6     	  # g/mol
Al_mass  = 27.0     	  # g/mol
T0	     = 300.           # K
Kb	     = 0.008314462175 # kJ/mol/K
N_a      = 6.02214076e23  # avogadro

NCbuff   = 10             # nlist update steps

NCout    = 1            # output profiles

tend	 = 10000e-3            # final time (ps)
dt	     = 1e-3           # timestep (ps)
NCsave   = 20 			  # save atoms


eV = 1.602176634
Na = 6.02214076
eV_Na = eV * Na

# atoms and domain
##############################################################
from Lattice import create_fcc

a = 3.61
supercell = (4,4,4)
ase_atoms = create_fcc("Cu", a, *supercell)


num_atoms = len(ase_atoms.numbers)
print("Number of Atoms: ", num_atoms)

# set temperature
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
MaxwellBoltzmannDistribution(ase_atoms, temperature_K=T0)
Stationary(ase_atoms)
ZeroRotation(ase_atoms)

from ase2mdcraft import convert_domain, convert_atoms

domain = convert_domain(ase_atoms, center = False)

Bxmin = domain.xmin
Bymin = domain.ymin
Bzmin = domain.zmin
Bxmax = domain.xmax
Bymax = domain.ymax
Bzmax = domain.zmax

print(Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax)

Lx = Bxmax - Bxmin
Ly = Bymax - Bymin
Lz = Bzmax - Bzmin
volume = Lx * Ly * Lz

print(Lx, Ly, Lz)

domain.set_periodic(0, True)
domain.set_periodic(1, True)
domain.set_periodic(2, True)

atoms = convert_atoms(ase_atoms)

atoms.data["uid"] = np.arange(1, atoms.data.size + 1)

# potential taken from MAISE database
##############################################################
potential = Potential("./Cu_model.pb", domain)
atoms.data["rcut"] = potential.rcut
print("rcut = ", atoms.data["rcut"])
##############################################################
# setup thermostat
thermostat = Langevin(
	beta		= 0.5, # ps^-1 
	temperature = T0,  # K 
	time_step   = dt,  # ps
	heat_x	    = 1,   # x axis heating
	heat_y	    = 1,   # y axis heating
	heat_z	    = 1,   # z axis heating
	Ux	        = 0,   # flow velocity x
	Uy	        = 0,   # flow velocity y
	Uz	        = 0    # flow velocity z
)

solver_cold = Solver(potential, threads=threads)
solver      = Solver(potential, thermostat, threads)

# allows to quickly turn off the thermostat
#solver      = Solver(potential, threads=threads)

stepper = Stepper(solver, threads=threads)
# setup neighbors list
atoms.data["rns"] = atoms.data["rcut"]
nlist = VerletList(atoms, atoms, domain=domain, threads=threads)
nlist.update()

##############################################################

# init profile for analysis
import analysis

nbins = 1
profile = analysis.Profile(
	xmin = Bxmin, xmax = Bxmax, 
	Ly   = Bymax - Bymin, 
	Lz   = Bzmax - Bzmin, 
	nbins = nbins
)

# prepare output
if not os.path.exists("./log/"):
	os.makedirs("./log/")
if not os.path.exists("./atoms/"):
	os.makedirs("./atoms/")

logfile = open("./log/log-prepare.dat","w")
logfile.write(
	"%06s %15s %15s %15s %15s %15s %15s %15s %15s\n" % (
		"step", "time", "potentialEnergy", "totalEnergy",
		"kineticEnergy", "temperature", "virial", "pressure", "density"
	)
)

# prepare main cycle
ICbuff   = 0
ICout    = 0
ICsave   = 0
ICwait   = 0
ICskip   = 0

time	 = 0.0
step	 = 0

start = Time.time()

# Calculate initial forces and energy
solver_cold.prepare(atoms, atoms, nlist.get())
solver_cold.forces(atoms,  atoms, nlist.get())
#solver_cold.virials(atoms, atoms, nlist.get())

def logout(atoms):
	energy_potential, \
	energy_total, \
	kinetic_energy, \
	temperature, \
	virial, \
	pressure, \
	density = analysis.averageValues(atoms, volume)
					
	elapsed = Time.time() - start
	logfile.write("%06d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n" %(
		step, time, 
		energy_potential, 
		energy_total, 
		kinetic_energy, 
		temperature, 
		virial,
		pressure, 
		density
	))
	print("step: %06d, phys_time: %6.2f ps, sys_time: %6.2f s" % (step, time, elapsed))
	logfile.flush()

logout(atoms)

##############################################################
# start main cycle
while step < int(tend/dt):

 	# here and in solver.forces are the same atoms
 	# prepare for all the atoms is done here
 	stepper.make_step(atoms, nlist.get(), dt) # here forces stage is done!

 	time += dt
 	step += 1

 	ICbuff = step % NCbuff
 	ICout  = step % NCout
 	ICsave = step % NCsave

 	if ICout  == 0:
 		#solver.virials(atoms, atoms, nlist.get())
 		logout(atoms)

 	if ICsave == 0:
 		save_file = hdf.File("atoms/prep%06d.h5" % step, 'w')
 		save_file["atoms"] = atoms.data
 		save_file["domain"] = np.array([Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax])
 		save_file.close()

 	if ICbuff == 0:
 		nlist.update()

logfile.close()
