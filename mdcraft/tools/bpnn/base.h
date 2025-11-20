#pragma once

#include <cmath>
#include <functional>

#include <mdcraft/tools/bpnn/autodiff_utils.h>

namespace mdcraft {
namespace bpnn {

using namespace ::mdcraft::bpnn::autodiff;

using TypeFc = std::function<double(double)>;

using TypeG2 = std::function<double(double)>;
using TypeG3 = std::function<double(double)>;
using TypeG4 = std::function<double(double, double, double, double)>;
using TypeG5 = std::function<double(double, double, double)>;

using TypedG4 = std::function<double(double, double, double, double, double, double, double)>;
using TypedG5 = std::function<double(double, double, double, double, double, double)>;

using TypeG4s = std::function<double(double, double, double, double, double, double, double, double, double)>;

class IBPDescriptor {
	using vectori2D = std::vector<std::vector<int>>;
	using vectord2D = std::vector<std::vector<double>>;
public:
	double      Rc()           const;
	double      Rs()           const;
	unsigned    numRadial()    const;
	unsigned    typeRadial()   const;
	unsigned    numAngular()   const;
	unsigned    typeAngular()  const;

	const std::vector<unsigned>&  shiftRadial()  const;
	const std::vector<unsigned>&  shiftAngular() const;

	double cutoff_function(double Rij) const;
	// each G function is called by its order number with specific MD arguments
	double G2_function(std::size_t ord, double Rij) const;
	double G3_function(std::size_t ord, double Rij) const;
	double G4_function(std::size_t ord, double Rij, double Rik, double Rjk, double cos) const;
	double G5_function(std::size_t ord, double Rij, double Rik, double cos) const;

	double G2_function_d(std::size_t ord, double Rij) const;
	double G3_function_d(std::size_t ord, double Rij) const;
	double G4_function_d(std::size_t ord, double Rij, double Rik, double Rjk,
						 double cos, double dcos, double dRij, double dRik) const;
	double G5_function_d(std::size_t ord, double Rij, double Rik,
						 double cos, double dcos, double dRij, double dRik) const;

	double G4_function_stress(std::size_t ord, double Rij, double Rik, double Rjk,
						 double cos, double dcos, double dRij, double dRik, double R_ij, double R_ik) const;

	double radial_function(std::size_t ord, double Rij) const;
	double angular_function(std::size_t ord, double Rij, double Rik, double Rjk, double cos) const;

	double radial_function_d(std::size_t ord, double Rij) const;
	double angular_function_d(std::size_t ord, double Rij, double Rik, double Rjk,
							  double cos, double dcos, double dRij, double dRik) const;
	double angular_function_stress(std::size_t ord, double Rij, double Rik, double Rjk,
							  double cos, double dcos, double dRij, double dRik, double R_ij, double R_ik) const;

protected:
	double    Rc_          {0.0};
	double    Rs_          {0.0};
	unsigned  numRadial_     {0};
	unsigned  numAngular_    {0};
	unsigned  typeRadial_    {0};
	unsigned  typeAngular_   {0};

	vectord2D 									   dzetas_;
	vectord2D									  lambdas_;
	vectord2D									   kappas_;

	// here and after assume an easy version for which
	// Rc_00 == Rc_01 == Rc_11 == Rc from second table in maise 'model' file
	TypeFc					   						   fc_;

	std::vector<TypeG2>						 G2iv_, dG2iv_;
	std::vector<TypeG3>						 G3iv_, dG3iv_;
	std::vector<TypeG4>						 		 G4iv_;
	std::vector<TypeG5>						 		 G5iv_;

	std::vector<TypedG4>						    dG4iv_;
	std::vector<TypedG5>						    dG5iv_;

	std::vector<TypeG4s>						    dG4iv_s_;

	std::vector<unsigned>		              shiftRadial_;                                     
	std::vector<unsigned>		             shiftAngular_;

	// functional patterns to be specified later for specific parameters
	static constexpr auto pattern_fc_ = [] (double Rij, double Rc) {
		assert(Rij < Rc);
		return 0.5 * ( std::cos(M_PI * Rij / Rc) + 1.0 );
	};

	// each Gi is assumed to handle only the part without multiplication by fc
	// each dGi contains both parts

	static constexpr auto pattern_G2_ = [] (double Rij, double Rs, double eta) {
		return std::exp( - eta * (Rij - Rs)*(Rij - Rs) );
	};

	static constexpr auto pattern_dG2_ = [] (double Rij, double Rs, double eta, double Rc) {
		assert(Rij < Rc);
		auto fc = [&] (double r) {
			return 0.5 * ( std::cos(M_PI * r / Rc) + 1.0 );
		};
		auto dfc = [&] (double r) {
			return (-0.5)*std::sin(M_PI * Rij / Rc) * M_PI / Rc;
		};
		return std::exp( - eta * (Rij - Rs)*(Rij - Rs) )
			 * ( - 2.0*eta*(Rij - Rs) * fc(Rij) + dfc(Rij) );
	};

	static constexpr auto pattern_G3_ = [] (double Rij, double kappa) {
		return std::cos(kappa * Rij);
	};

	static constexpr auto pattern_dG3_ = [] (double Rij, double kappa, double Rc) {
		assert(Rij < Rc);
		return std::cos(kappa * Rij) * (-0.5)*std::sin(M_PI * Rij / Rc) * M_PI / Rc
			 - std::sin(kappa * Rij) * kappa * 0.5 * ( std::cos(M_PI * Rij / Rc) + 1.0 );
	};

	static constexpr auto pattern_G4_ = [] (
		double Rij, double Rik, double Rjk, double cos_theta_ijk,
		double dzeta, double lambda, double eta
	) {
		return std::pow(2.0, 1.0 - dzeta) * std::pow(1.0 + lambda * cos_theta_ijk, dzeta)
		                       			  * std::exp( - eta * (Rij*Rij + Rik*Rik + Rjk*Rjk) );
	};

	static constexpr auto pattern_dG4_ = [] (
		double Rij, double Rik, double Rjk,
	 	double cos, double dcos, double dRij, double dRik,
		double dzeta, double lambda, double eta, double Rc
	) {
		assert((Rij < Rc) && (Rik < Rc) && (Rjk < Rc));
		auto fc = [&] (double r) {
			return 0.5 * ( std::cos(M_PI * r / Rc) + 1.0 );
		};
		auto dfc = [&] (double r) {
			return (-0.5)*std::sin(M_PI * r / Rc) * M_PI / Rc;
		};
		auto exp = [&] (double r1, double r2, double r3) {
			return std::exp( - eta * (r1*r1 + r2*r2 + r3*r3) );
		};
		return  fc(Rjk)
	   		   *exp(Rij, Rik, Rjk)
	   		   *std::pow(2.0, 1.0 - dzeta)
	   		   *std::pow(1.0 + lambda*cos, dzeta - 1.0)
		   	   *(  ( dzeta*lambda*dcos
		   	      + (1.0 + lambda*cos)*(-2.0*eta)*(Rij*dRij + Rik*dRik) ) * fc(Rij)*fc(Rik)
		   	      + (1.0 + lambda*cos)*(dRij*dfc(Rij)*fc(Rik) + dRik*dfc(Rik)*fc(Rij))
		   	    );
	};

	static constexpr auto pattern_dG4_s_ = [] (
		double Rij, double Rik, double Rjk,
	 	double cos, double dcos, double dRij, double dRik,
	 	double R_ij, double R_ik,
		double dzeta, double lambda, double eta, double Rc
	) {
		assert((Rij < Rc) && (Rik < Rc) && (Rjk < Rc));
		auto fc = [&] (double r) {
			return 0.5 * ( std::cos(M_PI * r / Rc) + 1.0 );
		};
		auto dfc = [&] (double r) {
			return (-0.5)*std::sin(M_PI * r / Rc) * M_PI / Rc;
		};
		auto exp = [&] (double r1, double r2, double r3) {
			return std::exp( - eta * (r1*r1 + r2*r2 + r3*r3) );
		};

		auto R_kj = R_ij - R_ik; // R_jk = Rj - Rk

		return  exp(Rij, Rik, Rjk)
	   		   *std::pow(2.0, 1.0 - dzeta)
	   		   *std::pow(1.0 + lambda*cos, dzeta - 1.0)
		   	   *(
		   	   		fc(Rik)*fc(Rjk) * ((-2*eta) * fc(Rij) * (1.0 + lambda*cos) - lambda*dzeta*cos*fc(Rij)/Rij/Rij + (1 + lambda*cos)*dfc(Rij)/Rij) * R_ij*R_ij +
		   	   		fc(Rij)*fc(Rjk) * ((-2*eta) * fc(Rik) * (1.0 + lambda*cos) - lambda*dzeta*cos*fc(Rik)/Rik/Rik + (1 + lambda*cos)*dfc(Rik)/Rik) * R_ik*R_ik +
		   	   		fc(Rik)*fc(Rij) * ((-2*eta) * fc(Rjk) * (1.0 + lambda*cos)                                    + (1 + lambda*cos)*dfc(Rjk)/Rjk) * R_kj*R_kj +
		   	   		2.0*lambda*dzeta * fc(Rij)*fc(Rjk)*fc(Rik)/Rij/Rik * R_ij*R_ik
		   	   	);
	};

	static constexpr auto pattern_G5_ = [] (
		double Rij, double Rik, double cos_theta_ijk,
		double dzeta, double lambda, double eta
	) {
		return std::pow(2.0, 1.0 - dzeta) * std::pow(1.0 + lambda * cos_theta_ijk, dzeta)
		       							  * std::exp( - eta * (Rij*Rij + Rik*Rik) );
	};

	static constexpr auto pattern_dG5_ = [] (
		double Rij, double Rik,
	 	double cos, double dcos, double dRij, double dRik,
		double dzeta, double lambda, double eta, double Rc
	) {
		assert((Rij < Rc) && (Rik < Rc));
		auto fc = [&] (double r) {
			return 0.5 * ( std::cos(M_PI * r / Rc) + 1.0 );
		};
		auto dfc = [&] (double r) {
			return (-0.5)*std::sin(M_PI * r / Rc) * M_PI / Rc;
		};
		auto exp = [&] (double r1, double r2) {
			return std::exp( - eta * (r1*r1 + r2*r2) );
		};
		return  exp(Rij, Rik)
	   		   *std::pow(2.0, 1.0 - dzeta)
	   		   *std::pow(1.0 + lambda*cos, dzeta - 1.0)
		   	   *(  ( dzeta*lambda*dcos
		   	      + (1.0 + lambda*cos)*(-2.0*eta)*(Rij*dRij + Rik*dRik) ) *fc(Rij)*fc(Rik)
		   	      + (1.0 + lambda*cos)*(dRij*dfc(Rij)*fc(Rik) + dRik*dfc(Rik)*fc(Rij))
		   	    );
	};
};

class IBPNN {
public:
	unsigned short                nKinds() const;
	const int                    nLayers() const;
	Tens2<double>&         		 biases( int kind, int layer);
	Tens2<double>&         		 weights(int kind, int layer);
	const Tens2<double>&         biases( int kind, int layer) const;
	const Tens2<double>&         weights(int kind, int layer) const;

	const long int  			 neurons(std::size_t i) const;

protected:
	unsigned short                       nKinds_;

	std::vector<int>                   architec_;

	// for each type
	struct BnW {
		std::vector<Tens2<double>>       biases_;
		std::vector<Tens2<double>>      weights_;
	};
	std::vector<BnW>					   anns_;
};

} // bpnn
} // mdcraft