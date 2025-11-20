#include <mdcraft/solver/thermostat/compound.h>

namespace mdcraft {
namespace solver {
namespace thermostat {

Compound::Compound() {}

void Compound::add(Base& thermostat) {
	thermostats.push_back(&thermostat);
}

void Compound::apply_one(Atoms::iterator atom) {
	for (auto& thermostat : thermostats) {
		thermostat->apply_one(atom);
	}
}

} // thermostat
} // solver
} // mdcraft