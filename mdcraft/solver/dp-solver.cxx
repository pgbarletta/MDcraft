#ifdef mdcraft_ENABLE_DeePMD

#include <mdcraft/solver/dp-solver.h>

namespace mdcraft {
namespace solver {

using ::mdcraft::tools::Threads;
using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

DP_Solver::DP_Solver(
	PotentialMB&  potential,
	Thermostat& thermostat,
	Threads& threads
) : ManyBodySolver(
		potential,
	    thermostat,
	    threads 
    )
{
}

void DP_Solver::prepare_one(
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

# endif
