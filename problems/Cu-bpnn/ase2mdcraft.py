import numpy as np

from ase import Atoms as AtomsASE
from mdcraft.data import Atoms

from mdcraft.lattice import Domain


def convert_domain(ase_atoms : AtomsASE, center : bool = True) -> Domain:
	Bxmin, Bymin, Bzmin = 0, 0, 0
	Bxmax, Bymax, Bzmax = ase_atoms.cell.lengths() * 0.1 # [A] -> [nm]
	if center:
		Bxmin -= Bxmax/2
		Bymin -= Bymax/2
		Bzmin -= Bzmax/2
		Bxmax -= Bxmax/2
		Bymax -= Bymax/2
		Bzmax -= Bzmax/2
	domain = Domain(
		Bxmin, Bxmax,
		Bymin, Bymax,
		Bzmin, Bzmax
	)
	return domain

def convert_atoms(ase_atoms : AtomsASE) -> Atoms:
	"""
		Mass is in atomic units in ASE,

		Atomic kinds are countered starting
		from zero in order, corresponding to
		atomic numbers.
		E.g., ['S', 'Zn'] -> [16, 30] (ASE) -> [0, 1] (mdcraft)
	"""
	ase_masses = ase_atoms.get_masses()
	size = len(ase_masses)
	atoms = Atoms(size)

	atoms.data["m"]      = ase_masses

	ase_kinds = ase_atoms.get_atomic_numbers()
	map_to_counter = lambda a: np.unique(a, return_inverse=True)[1]
	atoms.data["kind"] = map_to_counter(ase_kinds)

	positions  = ase_atoms.get_positions()
	# [A] -> [nm]
	atoms.data["r"]["x"] =  positions[:, 0] * 0.1
	atoms.data["r"]["y"] =  positions[:, 1] * 0.1
	atoms.data["r"]["z"] =  positions[:, 2] * 0.1

	velocities = ase_atoms.get_velocities()
	# [ase units: sqrt(eV/a.m.u.)] -> [nm / ps]
	from ase.units import fs
	atoms.data["v"]["x"] = velocities[:, 0] * fs * 1e2
	atoms.data["v"]["y"] = velocities[:, 1] * fs * 1e2
	atoms.data["v"]["z"] = velocities[:, 2] * fs * 1e2

	# assert np.allclose(VCM, 0), "Center of Mass Velocity is not Zero!"
	return atoms


def iconvert_atoms(atoms : Atoms, cell, pbc=True) -> AtomsASE:

	def map_mass_to_atomic_number(masses):
		from ase.data import atomic_masses
		atomic_numbers = []
		for mass in masses:
			atomic_number = np.argmin(np.abs(atomic_masses - mass))
			atomic_numbers.append(atomic_number)
		return np.array(atomic_numbers)

	def map_atomic_number_to_symbol(atomic_numbers_list):
		from ase.data import atomic_numbers
		symbol_dict = {v: k for k, v in atomic_numbers.items()}
		symbols = [symbol_dict[atomic_number] for atomic_number in atomic_numbers_list]
		return symbols
	
	masses = atoms.data["m"]
	kinds = map_mass_to_atomic_number(masses)
	symbols = map_atomic_number_to_symbol(kinds)
	ase_atoms = AtomsASE(symbols, cell=cell, masses=masses)

	positions  = ase_atoms.get_positions()
	# [nm] -> [A]
	positions[:, 0] = atoms.data["r"]["x"] * 10
	positions[:, 1] = atoms.data["r"]["y"] * 10
	positions[:, 2] = atoms.data["r"]["z"] * 10
	ase_atoms.set_positions(positions)

	velocities  = ase_atoms.get_velocities()
	# [nm / ps] -> [ase units: sqrt(eV/a.m.u.)]
	from ase.units import fs
	velocities[:, 0] = atoms.data["v"]["x"] / (fs * 1e2)
	velocities[:, 1] = atoms.data["v"]["y"] / (fs * 1e2)
	velocities[:, 2] = atoms.data["v"]["z"] / (fs * 1e2)
	ase_atoms.set_velocities(velocities)

	return ase_atoms

if __name__ == "__main__":
	atoms = AtomsASE('Au', cell=[10, 10, 10], pbc=True)
	domain = convert_domain(atoms, center=False)
	assert domain.xmin ==  0
	assert domain.xmax == 10
	domain = convert_domain(atoms)
	assert domain.xmin == -5
	assert domain.xmax ==  5


	from Lattice import create_zincblende, create_fcc
	ase_atoms = create_zincblende(["S", "Zn"], 5.42, 3, 3, 3)
	cell = ase_atoms.get_cell()

	#ase_atoms = create_fcc("Cu", 3.5, 5, 5, 5)

	from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
	temperature = 300.0
	MaxwellBoltzmannDistribution(ase_atoms, temperature_K=temperature, force_temp=True)
	Stationary(atoms)
	ZeroRotation(atoms)

	atoms = convert_atoms(ase_atoms)
	# print(atoms.data["m"])
	# print(atoms.data["kind"])
	# print(atoms.data["r"])

	# print(type(atoms.data["v"]))

	squared_v = np.empty_like(atoms.data["v"])
	for field in squared_v.dtype.fields.keys():
		squared_v[field] = atoms.data["v"][field] ** 2

	# tempx = squared_v["x"] * atoms.data["m"]
	# tempy = squared_v["y"] * atoms.data["m"]
	# tempz = squared_v["z"] * atoms.data["m"]
	
	# print(tempx + tempy + tempz)

	result = 0.0
	for field in squared_v.dtype.fields.keys():
		result += squared_v[field] * atoms.data["m"] #* 0.5 * 12.044 / atoms.data.size

	amu = 1.66053906892
	erg_to_kelvin = 7.2429716666
	# [nm / ps]**2 * amu -> Kelvin:
	kinerg_to_kelvin = amu * erg_to_kelvin * 10

	print("Temperature from atoms:")
	print(result.sum() * 0.5 / atoms.data.size * kinerg_to_kelvin / 1.5)
	print("Exprected temperature:")
	print(temperature)

	ase_atoms = iconvert_atoms(atoms, cell)
	from ase.visualize import view
	#view(ase_atoms)