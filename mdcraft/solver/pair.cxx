#include <mdcraft/solver/pair.h>

namespace mdcraft::solver {

Pair::Pair(
	Boundary&   boundary,
	Thermostat& thermostat,
	Threads&    threads
) : boundary   (boundary),
    thermostat (thermostat), 
    pool       (threads)
{
	boundary.bind_solver(*this);
}

void Pair::prepare(
	Atoms&     atoms,
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		prepare_one(atomit, neibs, nlist[i]);
		atom.bc = 0x0; boundary.prepare_one(atomit, neibs, nlist[i]);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::prepare(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
) {
	decomp.exchange_start();
	prepare(decomp.locals(), decomp.locals(), nlist1);
	decomp.exchange_end();
	prepare(decomp.locals(), decomp.aliens(), nlist2);
}
#endif

void Pair::charge(
	Atoms&     atoms, 
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		charge_one(atomit, neibs, nlist[i]);
		if (atom.bc) boundary.charge_one(atomit, neibs, nlist[i]);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::charge(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
	decomp.exchange_start();
	charge(decomp.locals(), decomp.locals(), nlist1);
	decomp.exchange_end();
	charge(decomp.locals(), decomp.aliens(), nlist2);
}
#endif

void Pair::forces(
	Atoms&     atoms, 
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		force_one(atomit, neibs, nlist[i]);
		thermostat.apply_one(atomit);
		if (atom.bc) boundary.force_one(atomit, neibs, nlist[i]);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::forces(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
	decomp.exchange_start();
	forces(decomp.locals(), decomp.locals(), nlist1);
	decomp.exchange_end();
	forces(decomp.locals(), decomp.aliens(), nlist2);
}
#endif

void Pair::virials(
	Atoms&     atoms, 
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		virial_one(atomit, neibs, nlist[i]);
		if (atom.bc) boundary.virial_one(atomit, neibs, nlist[i]);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::virials(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
	decomp.exchange_start();
	virials(decomp.locals(), decomp.locals(), nlist1);
	decomp.exchange_end();
	virials(decomp.locals(), decomp.aliens(), nlist2);
}
#endif

} // namespace mdcraft::solver
