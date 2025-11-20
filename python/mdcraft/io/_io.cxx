#include <mdcraft/configuration.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef mdcraft_ENABLE_MPI
#include <mpi4py/mpi4py.h>
#include <mdcraft/tools/network.h>
#endif

#include <mdcraft/data/atom.h>
#include <mdcraft/io/vtk/variables.h>
#include <mdcraft/io/vtk/pvd_file.h>

namespace py = pybind11;

using namespace mdcraft::io;
using mdcraft::data::Atoms;

PYBIND11_MODULE(_mdcraft_io, m) {
	py::class_<PvdFile>(m, "PvdFile")
	.def(py::init([](std::string filename, std::string directory) {
		return new PvdFile(filename, directory);
	}),
		py::arg("filename"),
		py::arg("dir") = "output"
	)
#ifdef mdcraft_ENABLE_MPI
	.def(py::init([](std::string filename, std::string directory,
		py::object  comm_obj) {
		PvdFile* pvd;
		if (comm_obj.is_none()) {
			pvd = new PvdFile(filename, directory);
		}
		else {
			MPI_Comm comm = ((PyMPICommObject*) comm_obj.ptr())->ob_mpi;
			pvd = new PvdFile(comm, filename, directory);
		}
		return pvd;
	}),
		py::arg("filename"),
		py::arg("dir")  = "output",
		py::arg("comm") = py::none()
	)
#endif
	.def("save", [](PvdFile& pvd, Atoms& atoms, double timestep) {
		pvd.save(atoms, timestep);
	},
		py::arg("atoms"),
		py::arg("timestep")
	)
	.def_property("variables", nullptr,
		[](PvdFile& pvd, py::list names) {
			pvd.variables.clear();
			for (size_t i = 0; i < names.size(); ++i) {
				pvd.variables += names[i].cast<std::string>();
			}
		});
}