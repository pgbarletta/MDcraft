#pragma once
#include <vector>
#include <string>
#include <mdcraft/data/vector.h>
#include <mdcraft/solver/potential/base.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

class Pair : public Base {
public:
	Pair(double        Rcutoff,
		const double* Vpar,
		int           Nr1,
		const double  dr
	);
	Pair(std::string filename);

	// potential U(r)
	double value(const vector r);
	double value(const double* r2, std::size_t size);
	void   value(const double* r2, double* Uout, std::size_t size);
	// void   value(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// virial U'(r)/r
	matrix virial(const vector r);
	matrix virial(const vector* r, const double* r2, std::size_t size);
	void   virial(const vector* r, const double* r2, matrix* Vout, std::size_t size);
	void   virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// force 
	vector force(const vector r);
	vector force(const vector* r, const double* r2, std::size_t size);
	void   force(const vector* r, const double* r2, vector* Fout, std::size_t size);
	void   force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);

private:
	// functions interpolation
	void init_spline(const double* y, double* z, double H1, std::size_t N, std::size_t ib);

	std::vector<double> Vpar, dVpar;
	int Nr1, N;
	double dr, dr1;
};

} // potential
} // solver
} // mdcraft