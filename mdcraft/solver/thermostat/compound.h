#pragma once
#include <mdcraft/data/atom.h>
#include <mdcraft/solver/thermostat/base.h>

namespace mdcraft { 
namespace solver { 
namespace thermostat {

using ::mdcraft::data::Atoms;

class Compound : public Base {
public:
	Compound();

	void add(Base& thermostat);
	void apply_one(Atoms::iterator atom);

private:

	std::vector<Base*> thermostats;

};

} // thermostat
} // solver
} // mdcraft