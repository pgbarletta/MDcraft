#include <mdcraft/solver/thermostat/base.h>

namespace mdcraft { 
namespace solver { 
namespace thermostat {

Base::Base() {}

void Base::apply_one(Atoms::iterator atom) { return; }

} // thermostat
} // solver
} // mdcraft