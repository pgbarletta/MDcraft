#include <mdcraft/solver/ternary.h>

namespace mdcraft {
namespace solver {

Ternary::Ternary(
	BaseTernary& potential,
	Thermostat&  thermostat,
	Threads&     threads
) : potential  (potential),
	thermostat (thermostat), 
    pool       (threads)
{
	stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
		forces(atoms, neibs, nlist);
	});
}

void Ternary::prepare(
	Atoms&     atoms,
	Atoms&     neibs,
	NeibsList& nlist
) {
	potential.reset_natoms(atoms.size());
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		if (atom.tag == 0) return;
		atom.f  = vector::Zero();
		atom.V  = matrix::Zero();
		atom.Ep = 0.0;
		atom.rcut = potential.rcut();
		potential.prepare(atomit, neibs, nlist);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

void Ternary::forces(
	Atoms&     atoms, 
	Atoms&     neibs, 
	const NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		// ternary potential holds a domain instance and applies periodics internally 
		if (atom.tag != 0) potential.force(atomit, neibs, nlist);
		thermostat.apply_one(atomit);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

void Ternary::virials(
	Atoms&     atoms, 
	Atoms&     neibs, 
	const NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		if (atom.tag != 0) potential.virial(atomit, neibs, nlist);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

} // solver
} // mdcraft
