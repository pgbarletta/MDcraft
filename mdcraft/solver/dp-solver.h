#pragma once

#include <mdcraft/solver/mb-solver.h>

namespace mdcraft {
namespace solver {

class DP_Solver : public ManyBodySolver {
public:
	DP_Solver(
		PotentialMB&  potential,
		Thermostat& thermostat = thermostat::dummy_thermostat,
		Threads& threads    = ::mdcraft::tools::dummy_pool
	);

	void prepare_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) override;
};


} // solver
} // mdcraft
