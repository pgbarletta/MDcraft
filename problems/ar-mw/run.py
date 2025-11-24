import numpy as np
import h5py as hdf
import matplotlib.pyplot as plt
import math, sys, os, time as Time
from mdcraft.tools import Threads
threads = Threads()
# data 
from mdcraft.data import Atoms
# output
from mdcraft.neibs import VerletList
# geometry
from mdcraft.lattice import Domain
# solver
from mdcraft.solver.potential  import LJs as Potential
from mdcraft.solver.thermostat import Langevin, Compound as CompoundThermostat
from mdcraft.solver.boundary   import Periodic
from mdcraft.solver            import Single as Solver
from mdcraft.solver.stepper    import Verlet as Stepper

from MW import MovingWindow

# physical params
Ar_mass  = 40	 # mass of argon atoms [g/mol]
He_mass  = 4	 # mass of helium atoms [g/mol]
Ar_part  = 0.8   # % of argon atoms
He_part  = 1.0 - Ar_part # % of helium atoms
Tin      = 120.  # inflow temperature [K]
Tout     = 2*Tin # outflow temperature [K]
uin      = -0.4  # inflow velocity (initial) [km/s]
uout     = 0.0   # outflow velocity (initial) [km/s]
xfront   = -10   # desired position of the shock front [nm]
dt	     = 4e-3  # timestep [ps]
time     = 0.0   # current time [ps]
tend	 = 1000  # final time [ps]
# technical params
step    = 0        # current step
NCbuff  = 10       # nlist update steps
Kbuff   = 1.3      # neighbor list buffer size
NCMW    = 5*NCbuff # MW update (insert inflow atoms, remove outflow atoms)
NCout   = 100      # output profiles
NCsave  = 2000     # save atoms
NxMulti = 8        # multiply initial sample
nbins   = 128      # number of bins along x axis to calculate physical quantities
# list of quantities to plot during simulation
plotQuantities = [
	'concentration',
	'density',
	'velocity',
	'temperature',
	'pressure'
]
# interatomic potential params for argon
Rcutoff = 0.8125 # nm
potential = Potential(
	aVr	 = 1.03120074782442750, # kJ/mol
	rVr	 = 0.33841043857528683, 
	Rcutoff = Rcutoff
)
# read prepared atoms
infile = hdf.File("./atoms.h5", "r")
buffer_atoms = infile["atoms"][:]
buffer_atoms["rns"] = Kbuff*buffer_atoms["rcut"]
# setup masses for Ar and He atoms
rndindex = np.random.choice([0, 1], buffer_atoms.size, p=[Ar_part, He_part])
iAr = np.where(rndindex == 0)
iHe = np.where(rndindex == 1)
buffer_atoms["m"][iAr] = Ar_mass
buffer_atoms["m"][iHe] = He_mass
# scale velocities
buffer_atoms["v"]["x"][iHe] *= np.sqrt(Ar_mass/He_mass)
buffer_atoms["v"]["y"][iHe] *= np.sqrt(Ar_mass/He_mass)
buffer_atoms["v"]["z"][iHe] *= np.sqrt(Ar_mass/He_mass)
# buffer domain sizes
Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax = infile["domain"][:]
xminbuf = Bxmin
xmaxbuf = Bxmax
Lxbuf   = Bxmax - Bxmin
bufsize = buffer_atoms.size
# allocate atoms (copy of buffer atoms NxMulti times along x axis)
atoms = Atoms(bufsize*NxMulti)
# new domain sizes
Bxmax =  NxMulti*0.5*Lxbuf
Bxmin = -NxMulti*0.5*Lxbuf
Lx	=  Bxmax - Bxmin
Ly	=  Bymax - Bymin
Lz	=  Bzmax - Bzmin
volume = Lx*Ly*Lz
# setup atoms data
for i in range(NxMulti):
	atoms.data[i*bufsize:(i+1)*bufsize] = buffer_atoms[:]
	atoms.data["r"]["x"][i*bufsize:(i+1)*bufsize] += Bxmin - xminbuf + i*Lxbuf
# setup computational domain
domain = Domain(
	1.1*Bxmin, 1.1*Bxmax, 
	    Bymin,     Bymax, 
	    Bzmin,     Bzmax
)
domain.set_periodic(1, True) # y period
domain.set_periodic(2, True) # z period
periodic = Periodic(domain)
# setup initial velocity (inflow)
atoms.data["v"]["x"] += uin
# setup moving window
MW = MovingWindow(
	buffer  = buffer_atoms, 
	bxmin   = xminbuf, 
	bxmax   = xmaxbuf, 
	xin     = Bxmax, 
	xout    = Bxmin,
	xfront  = xfront, 
	uin     = uin,
	uout    = uout
)
thermostats = MW.thermostats(Tin = Tin, Tout = Tout, dt = dt) # default values
# solver with thermostat
solver = Solver(
	potential  = potential, 
	boundary   = periodic, 
	thermostat = thermostats,
	threads	   = threads
)
# solver without thermostat
solver_cold = Solver(
	potential  = potential, 
	boundary   = periodic, 
	threads	   = threads
)
# setup MD stepper
stepper = Stepper(
	solver  = solver,
	threads = threads
)
# setup neighbors list
nlist = VerletList(
	atoms     = atoms,
	neighbors = atoms,
	domain    = domain,
	threads   = threads
)
nlist.update()
# init construction of profiles
import analysis

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
if not os.path.exists("./profiles/"):
	for quantity in plotQuantities:
		os.makedirs("./profiles/" + quantity)

logfile = open("./log/log-run.dat","w")
logfile.write(
	"%06s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n" % (
		"step", "natoms", "time", "xfront", "vshock", "vin", "vout", "potentialEnergy", "totalEnergy",
		"kineticEnergy", "temperature", "virial", "pressure", "density"
	)
)

# prepare main cycle
ICbuff = 0
ICMW   = 0
ICout  = 0
ICsave = 0
ICwait = 0
ICskip = 0

time   = 0.0
step   = 0

start = Time.time()

# Calculate initial forces and energy 
solver_cold.prepare(atoms,  atoms, nlist.get())
solver_cold.forces(atoms,  atoms, nlist.get())
solver_cold.virials(atoms, atoms, nlist.get())

# function to evaluate shock front position
def get_xfront(atoms):
	velocity = profile.analyze(atoms, smooth = 1.5)['velocity']
	xstep = (Bxmax - Bxmin)/velocity.size
	vaverage = np.mean(velocity)
	iout = np.where(velocity > vaverage)
	iin  = np.where(velocity < vaverage)
	vout = np.mean(velocity[iout])
	vin  = np.mean(velocity[iin])
	imid = np.where((velocity > vin + 0.2*(vout - vin))*(velocity < vout - 0.2*(vout - vin)))
	return np.mean(Bxmin + 0.5*xstep + imid[0]*xstep)

# start main cycle
while time < tend:

	MW.frozen_save(atoms)
	stepper.make_step(atoms, nlist.get(), dt)
	MW.frozen_step(atoms, dt)

	time += dt
	step += 1

	ICbuff = step % NCbuff
	ICMW   = step % NCMW
	ICout  = step % NCout
	ICsave = step % NCsave

	if ICout  == 0:
		solver.virials(atoms, atoms, nlist.get())
		energy_potential, \
		energy_total, \
		kinetic_energy, \
		temperature, \
		virial, \
		pressure, \
		density = analysis.averageValues(atoms, volume)

		elapsed = Time.time() - start
		logfile.write("%06d %15d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n" %(
			step, atoms.data.size, time, 
			MW.xfront,
			MW.vshock,
			MW.uin,
			MW.uout,
			energy_potential, 
			energy_total, 
			kinetic_energy, 
			temperature, 
			virial, 
			pressure, 
			density
		))
		print("step: %06d, time: %6.2f ps, elapsed: %6.2f s" % (step, time, elapsed))
		logfile.flush()
		sys.stdout.flush()

		# analysis
		data = profile.analyze(atoms, smooth = 1.5)
		for quantity in plotQuantities:
			plt.plot(profile.bins, data[quantity])
			plt.savefig("profiles/" + quantity + "/%06d.png" % step)
			plt.close()

	if ICMW == 0:
		MW.update(atoms, dt = NCMW*dt, xfront = get_xfront(atoms))

	if ICsave == 0:
		save_file = hdf.File("atoms/run%06d.h5" % step, 'w')
		save_file["atoms"] = atoms.data
		save_file["physics"] = np.array([
			dt,  time, tend,
			Tin, Tout, xfront,
			MW.bpos,
			MW.uin,
			MW.uout
		], dtype=np.float64)
		save_file["cycles"] = np.array([
			step,
			NCbuff,
			NCMW,
			NCout,
			NCsave
		], dtype=np.int64)
		save_file["domain"] = np.array([
			Bxmin, Bxmax, 
			Bymin, Bymax, 
			Bzmin, Bzmax
		])
		save_file.close()
		
	if ICbuff == 0:
		nlist.update()
