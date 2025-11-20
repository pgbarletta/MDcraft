#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <mpi4py/mpi4py.h>

#include <mdcraft/decomp/VD3.h>

struct Center3D {
    double x;
    double y;
    double z;
};

namespace py = pybind11;

using mdcraft::lattice::Domain;
using mdcraft::data::Atoms;
using mdcraft::tools::Threads;

using mdcraft::decomp::Decomp;
using mdcraft::decomp::VD3;

PYBIND11_MODULE(_mdcraft_decomp, m) {

PYBIND11_NUMPY_DTYPE(Center3D,
    x,
    y,
    z
);

py::class_<Decomp>(m, "Decomp")
	.def_property_readonly("rank", [](Decomp& D) -> int {
		return D.rank();
	})
    .def_property_readonly("size", [](Decomp& D) -> int {
        return D.size();
	})
	.def_property_readonly("locals", [](Decomp& D) -> Atoms& {
		return D.locals();
	}, py::return_value_policy::reference_internal)
	.def_property_readonly("aliens", [](Decomp& D) -> Atoms& {
		return D.aliens();
	}, py::return_value_policy::reference_internal)
	.def_property_readonly("border", [](Decomp& D) -> Atoms& {
		return D.border();
	}, py::return_value_policy::reference_internal)
	.def_property("measurer", nullptr, [](Decomp& D, const std::string& type) {
		D.set_measurer(type);
	})
	.def("exchange", [](Decomp& D) { D.exchange(); })
	.def("exchange_start", [](Decomp& D) { D.exchange_start(); })
	.def("exchange_end",   [](Decomp& D) { D.exchange_end(); })
	.def("update", [](Decomp& D, bool verbose) {
		D.update(verbose);
	},
		py::arg("verbose") = false
	)
	.def("prebalancing",   [](Decomp& D, int n_iters, bool verbose) {
		D.prebalancing(n_iters, verbose);
	},
		py::arg("n_iters") = 15,
		py::arg("verbose") = false
	)
;

py::class_<VD3, Decomp>(m, "VD3")
    .def(py::init([](
            py::object             comm_obj,
            Atoms&                 atoms,
            Domain&                domain,
            int                    dimension,
            py::array_t<Center3D>& centers,
            Threads&               pool,
            double                 mobility,
            double                 centroidal,
            double                 growth_rate,
            const std::string&     measurer
    ) {
    	// extract communicator from comm_obj
    	MPI_Comm comm = ((PyMPICommObject*) comm_obj.ptr())->ob_mpi;
    	// fill coords if needed
	    auto cppcenters = std::vector<mdcraft::data::vector>(centers.size());
	    if (centers.size() > 0) {
	    	auto pycenters = centers.data(); 
	    	for (int i = 0; i < centers.size(); ++i) {
	    		cppcenters[i].x() = pycenters[i].x;
	    		cppcenters[i].y() = pycenters[i].y;
	    		cppcenters[i].z() = pycenters[i].z;
	    	}
	    }
    	VD3* vd3 = new VD3(
       		comm,
		    atoms,
			domain,
			pool,
			{
				.dimension   = dimension,
        		.mobility    = mobility,
        		.centroidal  = centroidal,
        		.growth_rate = growth_rate
			},
			cppcenters
        );
    	vd3->set_measurer(measurer);
    	return vd3;
    }), py::arg("comm"), 
        py::arg("atoms"), 
        py::arg("domain"),
        py::arg("dimension") = 1,
        py::arg("centers") = py::array_t<Center3D>(),
        py::arg("threads") = mdcraft::tools::dummy_pool,
        py::arg("mobility") = 0.2,
        py::arg("centroidal") = 0.25,
        py::arg("growth_rate") = 0.02,
        py::arg("measurer") = "time"
    )
    .def_property_readonly("centers", [](VD3& vd3) -> py::array {
    	py::array_t<Center3D> centers(vd3.size());
    	auto pycenters  = centers.mutable_data();
    	auto cppcenters = vd3.centers(); 
    	for (int i = 0; i < vd3.size(); ++i) {
    		pycenters[i].x = cppcenters[i].x();
    		pycenters[i].y = cppcenters[i].y();
    		pycenters[i].z = cppcenters[i].z();
    	}
        return centers;
    });

}
