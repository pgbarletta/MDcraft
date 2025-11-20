#pragma once

#include <vector>

#include <mdcraft/data/atom.h>
#include <mdcraft/lattice/domain.h>
#include <mdcraft/neibs/verlet-list.h>

namespace mdcraft {
namespace solver {

class Pair;

namespace boundary { 

using ::mdcraft::data::vector;
using ::mdcraft::data::Atoms;
using NeibsListOne = ::mdcraft::neibs::ListOne;
using Solver = ::mdcraft::solver::Pair;
using ::mdcraft::lattice::Domain;
using ::mdcraft::lattice::dummy_domain;

static Atoms dummy_atoms;
static NeibsListOne dummy_indices_list;

/** 
	\brief \~russian Базовый тип граничных условий для метода MD.
	    \~english Base class for MD boundary condition. 
	    \~
*/
class Base {
public:

	/** 
		\brief 
			\~russian Конструктор класса.
			\~english Class constructor.
			\~
	*/
	Base(Domain& domain = dummy_domain);
	/** 
		\brief 
			\~russian Виртуальный деструктор класса.
			\~english Virtual class constructor.
			\~
	*/
	virtual ~Base() {}

	virtual void bind_solver(Solver& solver);

	virtual void prepare_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);
	
	virtual void charge_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);

	virtual void force_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);

	virtual void virial_one(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);

	Domain& get_domain();

	virtual Atoms prepare_ghosts(
		Atoms::iterator atom,
		Atoms&          neibs,
		NeibsListOne&   nlist
	);

protected:

	Solver* solver;
	Domain& domain;
};

static Base dummy_boundary;

} // boundary
} // solver
} // mdcraft
