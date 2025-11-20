#include <mdcraft/solver/boundary/base.h>
#include "base.h"

namespace mdcraft::solver::boundary {

Base::Base(Domain& domain)
	: domain(domain) {}

void Base::bind_solver(Solver& S) {
	solver = &S;
}

void Base::prepare_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {}

void Base::charge_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {}

void Base::force_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {}

void Base::virial_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {}

Domain& Base::get_domain() {
	return domain;
}

Atoms Base::prepare_ghosts(
	Atoms::iterator atom,
	Atoms &neibs,
	NeibsListOne &nlist)
{
	throw std::runtime_error("Base::prepare_ghosts is not implemented");
	return Atoms(); // TODO: implement this method
}

} // namespace mdcraft::solver::boundary
