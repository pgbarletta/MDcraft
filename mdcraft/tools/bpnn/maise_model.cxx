#include <cassert>

#include <mdcraft/tools/helpers.h>
#include <mdcraft/tools/bpnn/maise_model.h>
#include <mdcraft/tools/parser/maise_model.h>

namespace mdcraft {
namespace bpnn {

using namespace std::placeholders;

#ifdef DEBUG
using ::mdcraft::tools::print;
using ::mdcraft::tools::print;
#endif

using ::mdcraft::parser::MaiseModelParser;

MaiseBPDescriptor::MaiseBPDescriptor(std::string_view filename)
 : parser(std::make_unique<MaiseModelParser>(filename))
{
	parser->readHead();
	parser->readBody();
	auto nspc = parser->getHead()->nspc;

	kappas_.resize(nspc); dzetas_.resize(nspc); lambdas_.resize(nspc);

	for (std::size_t i = 0; i < nspc; ++i) {
		kappas_[i]  = parser->getBody()->table_vals[i][4];
		dzetas_[i]  = parser->getBody()->table_vals[i][5];
		lambdas_[i] = parser->getBody()->table_vals[i][6];
	}

	auto lin_factor = parser->getBody()->b2a[0]
					* parser->getBody()->scaling[0] // TODO: scaling not imposed here, but in a python script!
					* 0.1; // [A] -> [nm]
	auto inv_sq_factor = 1.0 / lin_factor / lin_factor;
	// params linear in distance
	Rc_ = parser->getBody()->table_vals[0][0][0] * lin_factor;
	Rs_ = parser->getBody()->table_vals[0][1][0] * lin_factor;

	fc_ = std::bind(pattern_fc_, _1, Rc_);

	// params inverse square in distance
	auto eta1_scaled = parser->getBody()->table_vals[0][2];
	for (auto& e : eta1_scaled) e *= inv_sq_factor;
	auto eta2_scaled = parser->getBody()->table_vals[0][3];
	for (auto& e : eta2_scaled) e *= inv_sq_factor;

	auto mapping = parser->getBody()->table_maps[0];
	auto num_descrs = mapping.size();
	assert(num_descrs > 0 && mapping[0].size() == 9);
	for (std::size_t i = 0; i < num_descrs; ++i) {
		auto order = mapping[i][1];
		switch (order) { // assumed that one and only one type of radial and angular funcs allowed
			case 2:
			{
				auto ieta = mapping[i][4];
				auto eta  = eta1_scaled[ieta];
				G2iv_.push_back(std::bind(pattern_G2_, _1, Rs_, eta));
				dG2iv_.push_back(std::bind(pattern_dG2_, _1, Rs_, eta, Rc_));
				typeRadial_ = 2;
				++numRadial_;
				break;
			}
			case 3:
			{
				auto ikappa = mapping[i][6];
				auto kappa  = kappas_[0][ikappa];
				G3iv_.push_back(std::bind(pattern_G3_, _1, kappa));
				dG3iv_.push_back(std::bind(pattern_dG3_, _1, kappa, Rc_));
				typeRadial_ = 3;
				++numRadial_;
				break;
			}
			case 4:
			{
				auto ieta = mapping[i][5];
				auto eta  = eta2_scaled[ieta];
				auto idz  = mapping[i][7];
				auto dz   = dzetas_[0][idz];
				auto il   = mapping[i][8];
				auto l    = lambdas_[0][il];
				G4iv_.push_back(std::bind(pattern_G4_, _1, _2, _3, _4, dz, l, eta));
				dG4iv_.push_back(std::bind(pattern_dG4_, _1, _2, _3, _4, _5, _6, _7, dz, l, eta, Rc_));
				dG4iv_s_.push_back(std::bind(pattern_dG4_s_, _1, _2, _3, _4, _5, _6, _7, _8, _9, dz, l, eta, Rc_));
				typeAngular_ = 4;
				++numAngular_;
				break;
			}
			case 5:
			{
				auto ieta = mapping[i][5];
				auto eta  = eta2_scaled[ieta];
				auto idz  = mapping[i][7];
				auto dz   = dzetas_[0][idz];
				auto il   = mapping[i][8];
				auto l    = lambdas_[0][il];
				G5iv_.push_back(std::bind(pattern_G5_, _1, _2, _3, dz, l, eta));
				dG5iv_.push_back(std::bind(pattern_dG5_, _1, _2, _3, _4, _5, _6, dz, l, eta, Rc_));
				typeAngular_ = 5;
				++numAngular_;
				break;
			}
			default:
				break;
		}
	}

	fillShiftTable();
}

void MaiseBPDescriptor::fillShiftTable() {
	shiftRadial_.resize(20, 0.0);
	shiftRadial_[0] = 0;
	shiftRadial_[1] = 1*numRadial_   + 1*numAngular_;
	shiftRadial_[2] = 2*numRadial_   + 3*numAngular_;
	shiftRadial_[3] = 3*numRadial_   + 6*numAngular_;

	shiftAngular_.resize(20, 0.0);
	shiftAngular_[ 0] = 1*numRadial_;
	shiftAngular_[ 1] = 2*numRadial_ + 1*numAngular_;
	shiftAngular_[ 2] = 2*numRadial_ + 2*numAngular_;
	shiftAngular_[ 4] = 3*numRadial_ + 3*numAngular_;
	shiftAngular_[ 8] = 3*numRadial_ + 4*numAngular_;
	shiftAngular_[ 5] = 3*numRadial_ + 5*numAngular_;
	shiftAngular_[ 9] = 4*numRadial_ + 6*numAngular_;
	shiftAngular_[18] = 4*numRadial_ + 7*numAngular_; 
	shiftAngular_[10] = 4*numRadial_ + 8*numAngular_;
	shiftAngular_[13] = 4*numRadial_ + 9*numAngular_; 
}

MaiseBPDescriptor::~MaiseBPDescriptor() {}

} // bpnn
} // mdcraft
