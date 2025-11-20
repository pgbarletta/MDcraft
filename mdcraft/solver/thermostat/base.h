#pragma once
#include <mdcraft/data/atom.h>

namespace mdcraft { 
namespace solver { 
namespace thermostat {

using ::mdcraft::data::Atoms;

const double Kb = 0.008314462175;// kJ/mol/K

class Base {
public:
	Base();
	double temperature();
	virtual ~Base() {}
	virtual void apply_one(Atoms::iterator atom);
};

static Base dummy_thermostat;

} // thermostat
} // solver
} // mdcraft