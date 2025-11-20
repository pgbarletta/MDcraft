#include <mdcraft/lattice/domain.h>
#include <mdcraft/solver/boundary/periodic.h>
// #include <mdcraft/solver/boundary/barrier.h>
// #include <mdcraft/solver/boundary/compound.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>

namespace py = pybind11;

using Domain   = ::mdcraft::lattice::Domain;
using Boundary = ::mdcraft::solver::boundary::Base;
using Periodic = ::mdcraft::solver::boundary::Periodic;
// using Barrier  = ::mdcraft::solver::boundary::Barrier;
// using Compound = ::mdcraft::solver::boundary::Compound;

PYBIND11_MODULE(_mdcraft_solver_boundary, m) {

py::class_<Boundary>(m, "Boundary");

py::class_<Periodic, Boundary>(m, "Periodic")
	.def(py::init<Domain&>());

// py::class_<Barrier, Boundary>(m, "Barrier")
// 	.def(py::init<double, double, double, double, double>(),
// 		py::arg("xpos"),
// 		py::arg("normal") = 1.0,
// 		py::arg("Ubndr") = 0.0,
// 		py::arg("Ubnda") = 0.0,
// 		py::arg("R1cut") = 0.0
// 	)
// 	.def("update_pos", [](Barrier& barrier, double xpos) {
// 		barrier.update_pos(xpos);
// 	})
// 	;

// py::class_<Compound, Boundary>(m, "Compound")
// 	.def(py::init<>())
// 	.def("add_boundary", [](Compound& compound, Boundary& bc) {
// 		compound.add_boundary(bc);
// 	})
// 	;

}