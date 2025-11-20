import glob

import numpy as np
import h5py as hdf
from scipy.ndimage import gaussian_filter1d

from mdcraft.data import Atoms
from ase2mdcraft import iconvert_atoms

from ase.io import read, write
from ase.visualize import view


Kb       = 0.008314462175# kJ/mol/K   R/1000
m_a      = 1.6605390666e-24  # g
N_a      = 6.02214076e23

def averageValues(atoms, volume):
	volume_nm3 = volume
	volume *= 1e-21
	vMCx = np.mean(atoms.data["v"]["x"])
	vMCy = np.mean(atoms.data["v"]["y"])
	vMCz = np.mean(atoms.data["v"]["z"])
	V2 = (atoms.data["v"]["x"]-vMCx)**2 + \
		 (atoms.data["v"]["y"]-vMCy)**2 + \
		 (atoms.data["v"]["z"]-vMCz)**2
	energy_potential = np.sum(atoms.data["Ep"])
	try:
		energy_potential += np.sum(atoms.data["Em"])
	except ValueError:
		pass
	two_kinetic_energies = np.sum(V2*atoms.data["m"])
	energy_total = two_kinetic_energies/2 + energy_potential
	temperature = two_kinetic_energies/Kb/3/atoms.data.size
	virial = np.sum(atoms.data["V"]["xx"] + \
		          atoms.data["V"]["yy"] + \
		          atoms.data["V"]["zz"])		
	pressure = (two_kinetic_energies - virial)/3.0/volume/N_a
	density = m_a*np.sum(atoms.data["m"])/volume
	return energy_potential, \
	       energy_total, \
	       two_kinetic_energies/2, \
		   temperature, \
		   virial, \
		   pressure, \
		   density

class Profile:

	def __init__(self, xmin, xmax, nbins, Ly, Lz):
		self.nbins = nbins
		self.bin_edges = np.linspace(xmin, xmax, nbins + 1)
		self.dx   = self.bin_edges[1] - self.bin_edges[0]
		self.bins = self.bin_edges[:-1] + 0.5*self.dx
		self.Vbin = self.dx*Ly*Lz*1e-21 # cm3

	def analyze(self, atoms, smooth = -1):
		data = np.zeros(self.nbins, dtype=[
			('concentration', 'f8'), # 1/cm3
			('density', 'f8'),       # g/cm3
			('temperature', 'f8'),   # K
			('pressure', 'f8'),      # GPa
			('velocity', 'f8'),       # km/s
			('energyint', 'f8'), # kJ/g
			('energypot', 'f8') # kJ/g

		])
		# evaluate concentration
		counts, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges
		)
		data['concentration'] = counts/self.Vbin
		# evaluate density
		mass, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = atoms.data["m"]
		)
		data['density'] = mass*m_a/self.Vbin
		# evaluate temperature and velocity
		binid = np.digitize(atoms.data["r"]["x"], bins = self.bin_edges, right=False)
		V2 = np.zeros(atoms.data.size, dtype=np.float64)
		for i in range(1, self.nbins + 1):
			idx = np.where(binid == i)
			# m calong x axis at bin i
			vMCx = 0.0
			if idx[0].size > 0:
				vMCx = np.sum(atoms.data["m"][idx]*atoms.data["v"]["x"][idx])
				vMCx /= np.sum(atoms.data["m"][idx])
			V2[idx] = (atoms.data["v"]["x"][idx] - vMCx)**2 + \
			           atoms.data["v"]["y"][idx]**2 + \
			           atoms.data["v"]["z"][idx]**2
			data['velocity'][i - 1] = vMCx

		twoEkin, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = V2*atoms.data["m"]
		)
		ige0 = np.where(counts > 0)
		data['temperature'][ige0] = twoEkin[ige0]/(counts[ige0] * Kb * 3)

		data['energypot'], tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = atoms.data["Ep"]
		)

		data['energyint'][ige0] = (data['energypot'][ige0] + twoEkin[ige0]/2.0)/mass[ige0]
		data['energypot'] /= mass[ige0]
		
		virials = atoms.data["V"]["xx"] + \
		          atoms.data["V"]["yy"] + \
		          atoms.data["V"]["zz"]
		V, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = virials
		)
		data['pressure'] = (twoEkin - V)/3.0/self.Vbin/N_a
		if smooth > 0:
			data['concentration'] = gaussian_filter1d(data['concentration'], sigma = smooth)
			data['density']       = gaussian_filter1d(data['density'],       sigma = smooth)
			data['temperature']   = gaussian_filter1d(data['temperature'],   sigma = smooth)
			data['pressure']      = gaussian_filter1d(data['pressure'],      sigma = smooth)
			data['velocity']      = gaussian_filter1d(data['velocity'],      sigma = smooth)
			data['energypot']     = gaussian_filter1d(data['energypot'],     sigma = smooth)
			data['energyint']     = gaussian_filter1d(data['energyint'],     sigma = smooth)
		return data

def extract_mdcraft_domain(hdf_file):
	with hdf.File(hdf_file, 'r') as load_file:
		domain_bounds = np.array(load_file["domain"])
		Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax = domain_bounds
	return (Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax)

def extract_mdcraft_atoms(hdf_file):
	with hdf.File(hdf_file, 'r') as load_file:
		atoms_data = np.array(load_file["atoms"])
		atoms = Atoms(atoms_data.size)
		atoms.data["m"] = atoms_data["m"]
		atoms.data["kind"] = atoms_data["kind"]
		atoms.data["r"] = atoms_data["r"]
		atoms.data["v"] = atoms_data["v"]
		atoms.data["V"] = atoms_data["V"]
		atoms.data["Ep"] = atoms_data["Ep"]
		atoms.data["Em"] = atoms_data["Em"]
	return atoms

def extract_ase_atoms(step, atoms_folder="./atoms"):		
	hdf_file = f"{atoms_folder}/prep{step:06d}.h5"
	atoms  = extract_mdcraft_atoms(hdf_file)
	domain = extract_mdcraft_domain(hdf_file)
	cell = [domain[1] - domain[0], domain[3] - domain[2], domain[5] - domain[4]]
	ase_atoms = iconvert_atoms(atoms, cell)
	return ase_atoms

def view_atomic_frame(step, atoms_folder="./atoms"):
	view(extract_ase_atoms(step, atoms_folder))

def create_movie(atoms_folder="./atoms", delta=1, format='traj', max_frame=np.iinfo(np.int32).max):
	frames = []
	
	step = 0
	for f in sorted(glob.glob(f"{atoms_folder}/prep*.h5")):
		if (step > max_frame):
			break
		step += delta
		atoms = extract_ase_atoms(step, atoms_folder)
		frames.append(atoms)

	print(f"Loaded {len(frames)} atomic frames from {atoms_folder} folder")
	write(f"atoms.{format}", frames)

def show_movie(movie="./atoms.traj"):
	frames = read(movie, index=":")
	view(frames)

def extract_density(hdf_file, mass, atoms_folder="./atoms", nbins=100, axis='x'):
	"""
		Provides density of atomic species with mass m
		along a gicen axis.
		Returns density in [g/cc] and concentration in [1/nm^3]
	"""
	if isinstance(hdf_file, int):
		hdf_file = f"{atoms_folder}/prep{hdf_file:06d}.h5"

	atoms = extract_mdcraft_atoms(hdf_file)
	domain = extract_mdcraft_domain(hdf_file)
	cell = [domain[1] - domain[0], domain[3] - domain[2], domain[5] - domain[4]]
	volume = cell[0] * cell[1] * cell[2]
	axis_to_index = {'x' : 0, 'y' : 1, 'z' : 2}


	bins = np.linspace(domain[2 * axis_to_index[axis]],
					   domain[2 * axis_to_index[axis] + 1], nbins + 1)
	binw = bins[1] - bins[0]
	density = np.zeros(nbins)
	concentration = np.zeros(nbins)

	ii = np.where(atoms.data["m"] == mass)

	masses    = atoms.data["m"][ii]
	positions = atoms.data["r"][ii]

	for i in range(nbins):
		i_bin = (positions[axis] >= bins[i]) & (positions[axis] < bins[i + 1])
		n_bin = len(masses[i_bin])
		m_bin = np.sum(masses[i_bin])
		V_bin = binw * volume / cell[axis_to_index[axis]]
		concentration[i] = n_bin / V_bin
		density[i] 		 = m_bin / V_bin

	bin_centers = 0.5 * (bins[:-1] + bins[1:])

	amu = 1.66053906892
	return bin_centers, density * 1e-3 * amu, concentration


def show_density_frame(hdf_file, plot_density=True, atoms_folder="./atoms", nbins=100, axis='x'):
	"""
		Plots instant density distribution along given axis in [g/cc] for each atom kind
		If plot_density=False, plots concentration in [1/nm^3]
	"""
	from ase2mdcraft import map_mass_to_atomic_symbol

	step = 0
	if isinstance(hdf_file, int):
		step = hdf_file
		hdf_file = f"{atoms_folder}/prep{hdf_file:06d}.h5"

	atoms    = extract_mdcraft_atoms(hdf_file)
	masses   = np.unique(atoms.data["m"])

	plt.figure(figsize=(8,6))
	for m in masses:
		bins, density, concentration = extract_density(step, m, atoms_folder, nbins, axis)
		symbol = map_mass_to_atomic_symbol([m])
		if plot_density:
			plt.plot(bins, density, label=f'{symbol}')
		else:
			plt.plot(bins, concentration, label=f'{symbol}')
	plt.xlabel("x axis, nm")
	if plot_density:
		plt.ylabel("density, g/cm$^3$")
	else:
		plt.ylabel("concentration, nm$^{-3}$")
	plt.legend()
	plt.show()


def create_density_gif(atoms_folder="./atoms", plot_density=True, nbins=100, axis='x'):
	from matplotlib.animation import FuncAnimation

	from ase2mdcraft import map_mass_to_atomic_symbol

	frames = sorted(glob.glob(f"{atoms_folder}/prep*.h5"))

	fig, ax = plt.subplots()
	ax.set_xlabel("x axis, nm")
	if plot_density:
		plt.ylabel("density, g/cm$^3$")
		ax.set_ylim(0, 15)
	else:
		plt.ylabel("concentration, nm$^{-3}$")
		ax.set_ylim(0, 1000)

	domain = extract_mdcraft_domain(frames[0])
	cell = [domain[1] - domain[0], domain[3] - domain[2], domain[5] - domain[4]]
	axis_to_index = {'x' : 0, 'y' : 1, 'z' : 2}
	ax.set_xlim(0, cell[axis_to_index[axis]])

	atoms  = extract_mdcraft_atoms(frames[0])
	masses = np.unique(atoms.data["m"])
	lines = {}
	for m in masses:
		symbol = map_mass_to_atomic_symbol([m])
		line, = ax.plot([], [], lw=2, label=f"{symbol}")
		lines[m] = line
	ax.legend()

	def update(frame):
		atoms  = extract_mdcraft_atoms(frame)
		masses = np.unique(atoms.data["m"])
		for m in masses:
			bins, density, concentration = \
				extract_density(frame, m, atoms_folder, nbins, axis)
			if plot_density:
				print('update!')
				lines[m].set_data(bins, density)
			else:
				lines[m].set_data(bins, concentration)
		return line,

	ani = FuncAnimation(fig, update, frames=frames, interval=30)#, blit=True)
	#ani.save("layer.bmp", writer="ffmpeg", fps=10)
	plt.show()

