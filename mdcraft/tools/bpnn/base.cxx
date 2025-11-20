#include <cassert>
#include <stdexcept>

#include <mdcraft/tools/bpnn/base.h>

namespace mdcraft {
namespace bpnn {

double IBPDescriptor::Rc() const {
	return Rc_;
}

double IBPDescriptor::Rs() const {
	return Rs_;
}

unsigned IBPDescriptor::typeRadial() const {
	return typeRadial_;
}

unsigned IBPDescriptor::typeAngular() const {
	return typeAngular_;
}

unsigned IBPDescriptor::numRadial() const {
	return numRadial_;
}

unsigned IBPDescriptor::numAngular() const {
	return numAngular_;
}

double IBPDescriptor::cutoff_function(double Rij) const {
	assert(fc_ && (Rij > 0.0));
	return fc_(Rij);
}

double IBPDescriptor::G2_function(std::size_t ord, double Rij) const {
	assert(Rij > 0.0);
	return G2iv_[ord](Rij);
}

double IBPDescriptor::G3_function(std::size_t ord, double Rij) const {
	assert(Rij > 0.0);
	return G3iv_[ord](Rij);
}

double IBPDescriptor::G4_function(
	std::size_t ord, double Rij, double Rik, double Rjk, double cos
) const {
	assert(Rij > 0.0); assert(Rik > 0.0); assert(Rjk > 0.0); assert(cos >= -1.0001); assert(cos <= 1.0001);
	return G4iv_[ord](Rij, Rik, Rjk, cos);
}

double IBPDescriptor::G5_function(
	std::size_t ord, double Rij, double Rik, double cos
) const {
	assert(Rij > 0.0); assert(Rik > 0.0); assert(cos >= -1.0001); assert(cos <= 1.0001);
	return G5iv_[ord](Rij, Rik, cos);
}

double IBPDescriptor::radial_function(std::size_t ord, double Rij) const {
	if (typeRadial_ == 2) {
		return G2_function(ord, Rij) * cutoff_function(Rij);
	} else if (typeRadial_ == 3) {
		return G3_function(ord, Rij) * cutoff_function(Rij);
	} else {
		throw std::runtime_error("Invalid number of radial BP function. Please use only <2> or <3> versions\n");
	}
}

double IBPDescriptor::G2_function_d(std::size_t ord, double Rij) const {
	assert(Rij > 0.0);
	return dG2iv_[ord](Rij);
}

double IBPDescriptor::G3_function_d(std::size_t ord, double Rij) const {
	assert(Rij > 0.0);
	return dG3iv_[ord](Rij);
}

double IBPDescriptor::G4_function_d(
	std::size_t ord, double Rij, double Rik, double Rjk,
	double cos, double dcos, double dRij, double dRik
) const {
	assert(Rij > 0); assert(Rik > 0); assert(Rjk > 0); assert(cos >= -1.0001); assert(cos <= 1.0001);
	return dG4iv_[ord](Rij, Rik, Rjk, cos, dcos, dRij, dRik);
}

double IBPDescriptor::G4_function_stress(
	std::size_t ord, double Rij, double Rik, double Rjk,
	double cos, double dcos, double dRij, double dRik, double R_ij, double R_ik
) const {
	assert(Rij > 0); assert(Rik > 0); assert(Rjk > 0); assert(cos >= -1.0001); assert(cos <= 1.0001);
	return dG4iv_s_[ord](Rij, Rik, Rjk, cos, dcos, dRij, dRik, R_ij, R_ik);
}

double IBPDescriptor::G5_function_d(
	std::size_t ord, double Rij, double Rik,
	double cos, double dcos, double dRij, double dRik
) const {
	assert(Rij > 0.0); assert(Rik > 0.0); assert(cos >= -1.0001); assert(cos <= 1.0001);
	return dG5iv_[ord](Rij, Rik, cos, dcos, dRij, dRik);
}

double IBPDescriptor::radial_function_d(std::size_t ord, double Rij) const {
	if (typeRadial_ == 2) {
		return G2_function_d(ord, Rij);
	} else if (typeRadial_ == 3) {
		return G3_function_d(ord, Rij);
	} else {
		throw std::runtime_error("Invalid number of radial BP function. Please use only <2> or <3> versions\n");
	}
}

double IBPDescriptor::angular_function(std::size_t ord, double Rij,
							double Rik, double Rjk, double cos
) const {
	if (typeAngular_ == 4) {
		return G4_function(ord, Rij, Rik, Rjk, cos) * 2.0 // why ?
		       * cutoff_function(Rij) * cutoff_function(Rik) * cutoff_function(Rjk);
	} else if (typeAngular_ == 5) {
		return G5_function(ord, Rij, Rik, cos)
			   * cutoff_function(Rij) * cutoff_function(Rik);
	} else {
		throw std::runtime_error("Invalid number of angular BP function. Please use only <4> or <5> versions\n");
	}
}

double IBPDescriptor::angular_function_d(std::size_t ord, double Rij,
		double Rik, double Rjk, double cos, double dcos, double dRij, double dRik
) const {
	if (typeAngular_ == 4) {
		return G4_function_d(ord, Rij, Rik, Rjk, cos, dcos, dRij, dRik);
	} else if (typeAngular_ == 5) {
		return G5_function_d(ord, Rij, Rik, cos, dcos, dRij, dRik);
	} else {
		throw std::runtime_error("Invalid number of angular BP function. Please use only <4> or <5> versions\n");
	}
}

double IBPDescriptor::angular_function_stress(std::size_t ord, double Rij,
		double Rik, double Rjk, double cos, double dcos, double dRij, double dRik, double R_ij, double R_ik
) const {
	if (typeAngular_ == 4) {
		return G4_function_stress(ord, Rij, Rik, Rjk, cos, dcos, dRij, dRik, R_ij, R_ik);
	}
	throw std::runtime_error("WHAT? #1529");
}

const std::vector<unsigned>& IBPDescriptor::shiftRadial()  const { return shiftRadial_;  }

const std::vector<unsigned>& IBPDescriptor::shiftAngular() const { return shiftAngular_; }

unsigned short IBPNN::nKinds() const {
	return nKinds_;
}

const long int IBPNN::neurons(std::size_t i) const {
	return architec_[i];
}

const int IBPNN::nLayers() const {
	return architec_.size();
}

const Tens2<double>& IBPNN::biases(int kind, int layer) const {
	return anns_[kind].biases_[layer];
}

const Tens2<double>& IBPNN::weights(int kind, int layer) const {
	return anns_[kind].weights_[layer];
}

Tens2<double>& IBPNN::biases(int kind, int layer) {
	return anns_[kind].biases_[layer];
}

Tens2<double>& IBPNN::weights(int kind, int layer) {
	return anns_[kind].weights_[layer];
}

} // bpnn
} // mdcraft