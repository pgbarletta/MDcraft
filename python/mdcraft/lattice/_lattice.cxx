#include <mdcraft/data/atom.h>
#include <mdcraft/lattice/domain.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using ::mdcraft::lattice::Domain;
using ::mdcraft::data::Atoms;

PYBIND11_MODULE(_mdcraft_lattice, m) {

py::class_<Domain>(m, "Domain")
    .def(py::init([](
            double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax
        ) {
        auto domain = new Domain(
            xmin, xmax,
            ymin, ymax,
            zmin, zmax
        );
        return domain;
    }),
        py::arg("xmin") = 0.0,
        py::arg("xmax") = 0.0,
        py::arg("ymin") = 0.0,
        py::arg("ymax") = 0.0,
        py::arg("zmin") = 0.0,
        py::arg("zmax") = 0.0,
        R"(Class constructor. Defines boundaries of the 
        computational domain.

        Parameters
        ----------

        xmin : double, optional
            minimal X in domain 
        xmax : double, optional
            maximal X in domain 
        ymin : double, optional
            minimal Y in domain 
        ymax : double, optional
            maximal Y in domain 
        zmin : double, optional
            minimal Z in domain 
        zmax : double, optional
            maximal Z in domain 

        Returns
        -------

        domain : Domain
            Domain class object.
        )"
    )
    .def("belong", [](Domain& domain, double x, double y, double z) -> bool {
        return domain.belong({x, y, z});
    },
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        R"(Determines whether the point is inside the 
        computational domain.

        Parameters
        ----------

        x : double, required
            x-coordinate. 
        y : double, required
            y-coordinate. 
        z : double, required
            z-coordinate.  

        Returns
        -------

        out : bool
            Whether the point is inside the computational
            domain.
        )"
    )
    .def("belong", [](Domain& domain, Atoms& atoms) -> py::array {
        auto belong = py::array_t<int>(atoms.size());
        auto cbelong = belong.mutable_data();
        for (std::size_t i = 0; i < atoms.size(); ++i)
            cbelong[i] = domain.belong(atoms[i].r);
        return belong;
    },
        py::arg("atoms"),
        R"(Determines whether the atoms are inside the 
        computational domain.

        Parameters
        ----------

        atoms : Atoms, required
            Atoms to be checked for belonging to the 
            computational domain.  

        Returns
        -------

        out : ndarray
            Boolean array which marks the atoms are inside
            the computational domain.
    )")
    .def("fit_in_period", [](Domain& domain,
        Atoms&  atoms
    ) -> void {
        for (auto& atom : atoms)
            domain.fit_in_period(atom.r);
    },
        py::arg("atoms"),
        R"(Adjusts the coordinates of the atoms to the 
        periodic boundary conditions.

        Parameters
        ----------

        atoms : Atoms, required
            Atoms to be Adjusted.  

    )")
    .def("nearest_distance", [](Domain& domain,
        Atoms& atoms,
        int    axis) -> py::array {
        auto distance = py::array_t<double>(atoms.size());
        auto cdistance = distance.mutable_data();
        for (std::size_t i = 0; i < atoms.size(); ++i)
            cdistance[i] = domain.nearest_distance(atoms[i].r, axis);
        return distance;
    },
        py::arg("atoms"),
        py::arg("axis"),
        R"(Determines nearest distances between atoms and domain boundaries.

        Parameters
        ----------

        atoms : Atoms, required

        Returns
        -------

        distance : ndarray
            Nearest distances to the domain boundaries.

    )")
    .def("set_periodic", [](Domain& domain, 
        int  axis,
        bool periodic
    ) -> void {
        domain.set_periodic(axis, periodic);
    },
        py::arg("axis"),
        py::arg("periodic") = true,
        R"(Sets periodic boundary conditions.

        Parameters
        ----------

        axis : int, required
            Number of axis to set periodic 
            boundary condition (x: 0, y: 1, z: 2).

        periodic : bool
            Turn on/off periodic condition
    )")
    .def("reshape", [](Domain& domain, 
        double xmin, double xmax,
        double ymin, double ymax,
        double zmin, double zmax
    ) -> void {
        domain.reshape(
            xmin, xmax,
            ymin, ymax,
            zmin, zmax
        );
    },
        py::arg("xmin") = 0.0,
        py::arg("xmax") = 0.0,
        py::arg("ymin") = 0.0,
        py::arg("ymax") = 0.0,
        py::arg("zmin") = 0.0,
        py::arg("zmax") = 0.0,
        R"(Redefines boundaries of the 
        computational domain.

        Parameters
        ----------

        xmin : double, optional
            minimal X in domain 
        xmax : double, optional
            maximal X in domain 
        ymin : double, optional
            minimal Y in domain 
        ymax : double, optional
            maximal Y in domain 
        zmin : double, optional
            minimal Z in domain 
        zmax : double, optional
            maximal Z in domain 
        )"
    )
    .def_property_readonly("xmin", [](Domain& domain) -> double {
        return domain.xmin();
    })
    .def_property_readonly("ymin", [](Domain& domain) -> double {
        return domain.ymin();
    })
    .def_property_readonly("zmin", [](Domain& domain) -> double {
        return domain.zmin();
    })
    .def_property_readonly("xmax", [](Domain& domain) -> double {
        return domain.xmax();
    })
    .def_property_readonly("ymax", [](Domain& domain) -> double {
        return domain.ymax();
    })
    .def_property_readonly("zmax", [](Domain& domain) -> double {
        return domain.zmax();
    })
    .def_property_readonly("xsize", [](Domain& domain) -> double {
        return domain.xsize();
    })
    .def_property_readonly("ysize", [](Domain& domain) -> double {
        return domain.ysize();
    })
    .def_property_readonly("zsize", [](Domain& domain) -> double {
        return domain.zsize();
    })
    .def_property_readonly("xperiodic", [](Domain& domain) -> bool {
        return domain.xperiodic();
    })
    .def_property_readonly("yperiodic", [](Domain& domain) -> bool {
        return domain.yperiodic();
    })
    .def_property_readonly("zperiodic", [](Domain& domain) -> bool {
        return domain.zperiodic();
    })
    ;
}