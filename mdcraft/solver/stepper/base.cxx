#include <mdcraft/solver/stepper/base.h>

namespace mdcraft::solver::stepper {

Base::Base(
    ISolver& solver,
    Threads& threads
) : solver(solver),
    pool(threads)
{}

} // namespace mdcraft::solver::stepper
