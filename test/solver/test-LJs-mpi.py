from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("load mpi4py: rank %d, size %d." % (rank, size))

import numpy as np
import sys, os, time as Time
import h5py as hdf
from mdcraft.tools import Threads
threads = Threads(comm=comm)

from mdcraft.io import PvdFile
# data 
from mdcraft.data import Atoms
# output
from mdcraft.neibs import VerletList
# geometry
from mdcraft.lattice import Domain
from mdcraft.decomp import VD3
# solver
from mdcraft.solver.potential  import LJs as Potential
from mdcraft.solver.thermostat import Langevin
from mdcraft.solver.boundary   import Periodic
from mdcraft.solver            import Single as Solver
from mdcraft.solver.stepper    import Verlet as Stepper

# times
dt      = 4e-3   # ps

tend    = 300.0  # ps
kBuf    = 1.25   # buffer size (default 1.25)
NCbuff  = 10     # nlist update
NCout   = 10     # analyze data
NCsave  = 300    # save atoms
# atoms
Rcutoff = 0.8125 # nm
potential = Potential(
	aVr     = 1.03120074782442750, # kJ/mol
	rVr     = 0.33841043857528683, 
	Rcutoff = Rcutoff
)

Ar_mass = 39.948 # Argon g/mol
T0      = 120    # K
Kb      = 0.008314462175# kJ/mol/K
N_a     = 6.02214076 * 10**23

# setup atoms mesh
box_sizes = None
sample_sizes = None

# create on master process
atoms = Atoms(0)

if rank == 0:
	from mesh import Mesh
	mesh = Mesh(
		Lx   = 24.0,     # nm
		Ly   = 24.0,     # nm
		Lz   = 24.0,     # nm
		D    = 0.40323, # nm, from potential file
		kind = 'fcc'
	)
	atoms = mesh.setup()

	atoms.data["rcut"] = potential.rcut
	atoms.data["rns"] = kBuf*atoms.data["rcut"]
	atoms.data["m"] = Ar_mass

	# initial velocity distribution
	mu = 0.0
	sigma = (1.2*Kb*T0/Ar_mass)**0.5
	Vx = np.random.normal(mu, sigma, atoms.data.size)
	Vy = np.random.normal(mu, sigma, atoms.data.size)
	Vz = np.random.normal(mu, sigma, atoms.data.size)
	V_comx = np.sum(Vx)/atoms.data.size
	V_comy = np.sum(Vy)/atoms.data.size
	V_comz = np.sum(Vz)/atoms.data.size
	atoms.data["v"]["x"] = Vx - V_comx
	atoms.data["v"]["y"] = Vy - V_comy
	atoms.data["v"]["z"] = Vz - V_comz

	# setup boundary conditions
	box_sizes = mesh.box_sizes()
	sample_sizes = mesh.sample_sizes()

# broadcast sizes
box_sizes = comm.bcast(box_sizes, root=0)
sample_sizes = comm.bcast(sample_sizes, root=0)

Lx, Ly, Lz = box_sizes
Bxmin, Bymin, Bzmin, \
Bxmax, Bymax, Bzmax = sample_sizes

# create domain
volume = Lx * Ly * Lz
domain = Domain(
	Bxmin, Bxmax,
	Bymin, Bymax,
	Bzmin, Bzmax
)
domain.set_periodic(0, True)
domain.set_periodic(1, True)
domain.set_periodic(2, True)
periodic = Periodic(domain)

# Сжимает атомы в центр области
def transform(atoms, p):
	X = 2.0 * (atoms.data["r"]["x"]- Bxmin) / (Bxmax - Bxmin) - 1.0
	Y = 2.0 * (atoms.data["r"]["y"]- Bymin) / (Bymax - Bymin) - 1.0
	Z = 2.0 * (atoms.data["r"]["z"]- Bzmin) / (Bzmax - Bzmin) - 1.0

	atoms.data["r"]["x"] = 0.5 * ((Bxmin + Bxmax) + (Bxmax - Bxmin) * np.sign(X) * np.abs(X)**p)
	atoms.data["r"]["y"] = 0.5 * ((Bymin + Bymax) + (Bymax - Bymin) * np.sign(Y) * np.abs(Y)**p)
	atoms.data["r"]["z"] = 0.5 * ((Bzmin + Bzmax) + (Bzmax - Bzmin) * np.sign(Z) * np.abs(Z)**p)

if rank == 0:
	transform(atoms, 1.02)

# create decomposition (auto redistribute/update)
VD = VD3(
	comm      = comm,     # required
	atoms     = atoms,    # required
	domain    = domain,   # required
	dimension = 2,        # optional
	threads   = threads,  # optional
	centroidal= 0.25,     # optional
	measurer  = "size"    # optional
)
VD.prebalancing(verbose=True)
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
# solver with thermostat
solver = Solver(
	potential  = potential, 
	boundary   = periodic, 
	# thermostat = thermostat, # uncomment to use
	threads    = threads
)
# solver without thermostat
solver_cold = Solver(
	potential  = potential, 
	boundary   = periodic, 
	threads    = threads
)
# setup MD stepper
stepper = Stepper(
	solver  = solver,
	threads = threads
)
# setup neighbors list (locals, locals)
nlist1 = VerletList(
	atoms     = VD.locals,
	neighbors = atoms,
	domain    = domain,
	threads   = threads
)
nlist1.update()
# setup neighbors list (locals, aliens)
nlist2 = VerletList(
	atoms     = VD.locals,
	neighbors = VD.aliens,
	domain    = domain,
	threads   = threads
)
nlist2.update()

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
if not os.path.exists("./log/") and rank == 0:
	os.makedirs("./log/")
if not os.path.exists("./atoms/") and rank == 0:
	os.makedirs("./atoms/")

logfile = open("./log/log-prepare.pt%d.dat" % rank, "w")
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
solver_cold.prepare(VD, nlist1.get(), nlist2.get())
solver_cold.forces (VD, nlist1.get(), nlist2.get())
solver_cold.virials(VD, nlist1.get(), nlist2.get())

pvd = PvdFile("atoms", comm=comm)
pvd.variables = ["rank", "Ep", "rns", "rcut", "bc", "f"]
pvd.save(VD.locals, 0.0)

# start main cycle
while step < NCsave:

	stepper.make_step(VD, nlist1.get(), nlist2.get(), dt)

	time += dt
	step += 1

	ICbuff = step % NCbuff
	ICout  = step % NCout
	ICsave = step % NCsave

	if ICout  == 0:
		solver.virials(VD, nlist1.get(), nlist2.get())
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
		if rank == 0:
			print("step: %06d, phys_time: %6.2f ps, sys_time: %6.2f s" % (step, time, elapsed))
		logfile.flush()
		pvd.save(VD.locals, time)

	if ICsave == 0:
		save_file = hdf.File("atoms/prep%06d.pt%d.h5" % (step, rank), 'w')
		save_file["atoms"] = atoms.data
		save_file["domain"] = np.array([Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax])
		save_file.close()
		
	if ICbuff == 0:
		VD.update(verbose=False)
		nlist1.update()
		nlist2.update()

logfile.close()
