#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>

#include <mdcraft/solver/thermostat/langevin.h>
#include <mdcraft/solver/thermostat/compound.h>

namespace py = pybind11;

using Thermostat = ::mdcraft::solver::thermostat::Base;
using Langevin   = ::mdcraft::solver::thermostat::Langevin;
using Compound   = ::mdcraft::solver::thermostat::Compound;

PYBIND11_MODULE(_mdcraft_solver_thermostat, m) {

py::class_<Thermostat>(m, "Thermostat")
	.def(py::init<>());

py::class_<Langevin, Thermostat>(m, "Langevin")
	.def(py::init<double, double, double, int, int, int, double, double, double, double, double>(),
		py::arg("beta"),
		py::arg("temperature"),
		py::arg("time_step"),
		py::arg("heat_x") =  1,
		py::arg("heat_y") =  1,
		py::arg("heat_z") =  1,
		py::arg("Ux")     =  0.0,
		py::arg("Uy")     =  0.0,
		py::arg("Uz")     =  0.0,
		py::arg("xmin")   = -1e+15,
		py::arg("xmax")   =  1e+15
	)
	.def_property_readonly("temperature", [](Langevin& thermostat) -> double {
		return thermostat.temperature();
	})
	.def_property_readonly("friction", [](Langevin& thermostat) -> double {
		return thermostat.friction();
	})
	.def_property_readonly("Ux", [](Langevin& thermostat) -> double {
		return thermostat.average_velocity().x();
	})
	.def_property_readonly("Uy", [](Langevin& thermostat) -> double {
		return thermostat.average_velocity().y();
	})
	.def_property_readonly("Uz", [](Langevin& thermostat) -> double {
		return thermostat.average_velocity().z();
	})
	.def_property_readonly("xmin", [](Langevin& thermostat) -> double {
		return thermostat.xmin();
	})
	.def_property_readonly("xmax", [](Langevin& thermostat) -> double {
		return thermostat.xmax();
	})
	.def("set_average_velocity", [](Langevin& thermostat, 
		double Ux,  
		double Uy,
		double Uz
	) -> void {
		thermostat.set_average_velocity({Ux, Uy, Uz});
	},
		py::arg("Ux") =  0.0,
		py::arg("Uy") =  0.0,
		py::arg("Uz") =  0.0
	)
	.def("set_temperature", [](Langevin& thermostat, 
		double T
	) -> void {
		thermostat.set_temperature(T);
	})
	.def("set_friction", [](Langevin& thermostat, 
		double beta
	) -> void {
		thermostat.set_friction(beta);
	})
	.def("set_heating_axes", [](Langevin& thermostat, 
		int heat_x,
		int heat_y,
		int heat_z
	) -> void {
		thermostat.set_heating_axes(
			heat_x,
			heat_y,
			heat_z
		);
	},
		py::arg("heat_x") = 1,
		py::arg("heat_y") = 1,
		py::arg("heat_z") = 1
	)
	.def("set_dimensions", [](Langevin& thermostat, 
		double xmin,
		double xmax
	) -> void {
		thermostat.set_dimensions(
			xmin,
			xmax
		);
	},
		py::arg("xmin") = -1e+15,
		py::arg("xmax") =  1e+15
	)
	;

py::class_<Compound, Thermostat>(m, "Compound")
	.def(py::init<>())
	.def("add", [](Compound& compound, Thermostat& thermostat) -> void {
		compound.add(thermostat);
	});
}
