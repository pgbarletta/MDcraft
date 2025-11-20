#pragma once

#include <string>

#include <mdcraft/configuration.h>

#ifdef mdcraft_ENABLE_DeePMD

#include <mdcraft/solver/potential/base.h>

#include "deepmd/deepmd.hpp"

namespace mdcraft {
namespace solver {
namespace potential {

namespace deepmd_compat = deepmd::hpp;

class DeepModelingPotential : public BaseManyBody {
public:
	DeepModelingPotential(const std::string& filename, Domain& domain);
	~DeepModelingPotential();

	//void prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) override;
	void force(   Atoms& atoms, Atoms& neibs, NeibsList& nlist) override;
	// void virial(  Atoms& atoms, Atoms& neibs, NeibsList& nlist) override;

private:
	deepmd_compat::DeepPot deep_pot;

	unsigned numb_models;
	double cutoff;
	int numb_types;
	int numb_types_spin;
	int dim_fparam;
	int dim_aparam;

	std::vector<double> fparam;
	std::vector<double> aparam;
};

} // potential
} // solver
} // mdcraft

#endif