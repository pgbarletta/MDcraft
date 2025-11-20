#pragma once

#include <string_view>

#include <mdcraft/configuration.h>

#ifdef mdcraft_ENABLE_BPNN

#include <mdcraft/tools/bpnn/maise_ann.h>
#include <mdcraft/tools/bpnn/maise_model.h>
#include <mdcraft/tools/bpnn/autodiff_utils.h>

#include <mdcraft/tools/threads.h>

#include <mdcraft/solver/potential/base.h>

namespace mdcraft {
namespace solver {
namespace potential {

using ::mdcraft::bpnn::MaiseBPNN;
using ::mdcraft::bpnn::MaiseBPDescriptor;
using namespace ::mdcraft::bpnn::autodiff;

class MaiseBPNNPotential : public BaseTernary {
public:
	MaiseBPNNPotential(std::string_view filename, Domain& domain);
	~MaiseBPNNPotential();

	void prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) override;
	void force(   Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) override;
	void virial(  Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) override;

	void reset_natoms(std::size_t natoms) override;

private:
	MaiseBPNN                                                      ann_;
	MaiseBPDescriptor                                       descriptor_;

	static const long int INPUT_LAYER = 300;
	static const long int NEIBS_ESTIM = 512;

	struct GradsThreadInfo {
		Eigen::TensorFixedSize<double, Eigen::Sizes<INPUT_LAYER, 1>>     			  g;
		Eigen::TensorFixedSize<double, Eigen::Sizes<INPUT_LAYER, 1>>  			   dedg;
		Eigen::TensorFixedSize<double, Eigen::Sizes<INPUT_LAYER, 3, NEIBS_ESTIM>>  dgdr;
		Eigen::TensorFixedSize<double, Eigen::Sizes<INPUT_LAYER, 3, NEIBS_ESTIM>>  dgdr_X_r;
	};

	std::vector<GradsThreadInfo>								 tinfo_;

	// cite energy and gradients
	void dEdG(Atoms::iterator atomit, Atoms& atoms);

	// utility functions for creating input layer info and inner gradients (dG/dr)
	void refresh(Atoms::iterator atomit, Atoms& atoms, const NeibsList& nlist);
	void accumulate_inner_grads(Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist);
};

} // potential
} // solver
} // mdcraft

#endif