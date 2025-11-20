#include <mdcraft/solver/mb-solver.h>

namespace mdcraft {
namespace solver {

using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

ManyBodySolver::ManyBodySolver(
	PotentialMB&  potential,
	Thermostat& thermostat,
	Threads& threads
) : potential  (potential),
    thermostat (thermostat),
    pool       (threads)
{
	stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
		forces(atoms, neibs, nlist);
	});
}

void ManyBodySolver::prepare(
	Atoms&     atoms,
	Atoms&     neibs,
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		prepare_one(atomit, neibs, nlist[i]);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

void ManyBodySolver::forces(
	Atoms&     atoms,
	Atoms&     neibs,
	NeibsList& nlist
) {
	potential.force(atoms, neibs, nlist);

	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		thermostat.apply_one(atomit);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

void ManyBodySolver::prepare_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	atom.f  = vector::Zero();
	atom.V  = matrix::Zero();
	atom.Ep = 0.0;
	atom.rcut = potential.rcut();
}

} // solver
} // mdcraft
