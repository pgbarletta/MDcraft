#include <mdcraft/solver/isolver.h>


namespace mdcraft::solver {

ISolver::~ISolver() = default;

void ISolver::prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) { }

#ifdef mdcraft_ENABLE_MPI
void ISolver::prepare(Decomp& decomp, NeibsList& nlist1, NeibsList& nlist2) { }
#endif
} // namespace mdcraft::solver