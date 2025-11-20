#pragma once

#include <mdcraft/solver/isolver.h>

#include <mdcraft/tools/threads.h>
#include <mdcraft/solver/potential/base.h>
#include <mdcraft/solver/thermostat/base.h>

namespace mdcraft {
namespace solver {

using PotentialMB  = potential::BaseManyBody;
using Thermostat = thermostat::Base;
using tools::Threads;

using NeibsListOne = ::mdcraft::neibs::ListOne;

class ManyBodySolver : public ISolver {
public:
	ManyBodySolver(
		PotentialMB&  potential,
		Thermostat& thermostat = thermostat::dummy_thermostat,
		Threads& threads = tools::dummy_pool
	);

	virtual ~ManyBodySolver() {};

	void prepare(
		Atoms&     atoms,
		Atoms&     neibs,
		NeibsList& nlist
	) override;

	void forces(
		Atoms&     atoms,
		Atoms&     neibs,
		NeibsList& nlist
	);

	virtual void prepare_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);

protected:
	PotentialMB&  potential;
	Thermostat& thermostat;
	Threads& pool;
};

} // solver
} // mdcraft
