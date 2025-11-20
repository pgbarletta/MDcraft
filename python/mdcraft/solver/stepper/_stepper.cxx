#include <mdcraft/configuration.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/data/atom.h>
#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/solver/stepper/base.h>
#include <mdcraft/solver/stepper/verlet.h>
#ifdef mdcraft_ENABLE_MPI
#include <mdcraft/decomp/decomp.h>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using mdcraft::tools::Threads;
using mdcraft::data::Atoms;
using NeibsList = mdcraft::neibs::List;
using Stepper = mdcraft::solver::stepper::Base;
using mdcraft::solver::ISolver;
using mdcraft::solver::stepper::Verlet;
#ifdef mdcraft_ENABLE_MPI
using mdcraft::decomp::Decomp;
#endif

PYBIND11_MODULE(_mdcraft_solver_stepper, m) {

py::class_<Stepper>(m, "Stepper")
	.def("make_step", [](
		Stepper&    stepper, 
		Atoms&      atoms, 
		NeibsList&  nlist, 
		double      dt
	) -> void {
		stepper.make_step(atoms, nlist, dt);
	},
		py::arg("atoms"),
		py::arg("nlist"),
		py::arg("dt")
	)
#ifdef mdcraft_ENABLE_MPI
	.def("make_step", [](
		Stepper&    stepper,
		Decomp&     decomp,
		NeibsList&  nlist1,
		NeibsList&  nlist2,
		double      dt
	) -> void {
		stepper.make_step(decomp, nlist1, nlist2, dt);
	},
		py::arg("decomp"),
		py::arg("nlist1"),
		py::arg("nlist2"),
		py::arg("dt")
	)
#endif
	;

py::class_<Verlet, Stepper>(m, "Verlet")
	.def(py::init([](
		ISolver& solver,
		Threads& pool
	) {
		auto stepper = new Verlet(
			solver,
			pool
		);
		return stepper;
	}), py::arg("solver"),
	    py::arg("threads") = mdcraft::tools::dummy_pool
	)
	;
}
