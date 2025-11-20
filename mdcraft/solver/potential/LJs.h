#pragma once
#include <mdcraft/solver/potential/base.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

// pair Lennard-Jones potential V(r)[kJ/mol]
// X =  (r/rVr)^2 - X2min
// V(r) = 4* aVr * ( [rVr/r]^12 - [rVr/r]^6 + [aLJ3 + bLJ2*X]* X*X )

class LJs : public Base {
public:
	LJs(
		double aVr,
		double rVr,
		double Rcutoff
	);
	// potential U(r)
	double value(const vector r);
	double value(const double* r2, std::size_t size);
	void   value(const double* r2, double* Uout, std::size_t size);
	// void   value(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);
	// virial U'(r)/r
	matrix virial(const vector r);
	matrix virial(const vector* r, const double* r2, std::size_t size);
	void   virial(const vector* r, const double* r2, matrix* Vout, std::size_t size);
	void   virial(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);
	// force 
	vector force(const vector r);
	vector force(const vector* r, const double* r2, std::size_t size);
	void   force(const vector* r, vector* Fout, std::size_t size);
	void   force(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);

private:
	
	double aVr, aVr2, rVr, rVr1, rVr2;
	double aLJ, aLJ3, bLJ, bLJ2;

	const double X2min = 1.2599210498948732;

};

} // potential
} // solver
} // mdcraft