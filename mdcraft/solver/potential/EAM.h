#pragma once
#include <vector>
#include <string>
#include <mdcraft/data/vector.h>
#include <mdcraft/solver/potential/base.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

class EAM : public Base {
public:
	EAM(double        Rcutoff,
		const double* Femb,
		const double* Rho,
		const double* Vpar,
		int           Nr1,
		const double  dr,
		int           Nrh1,
		const double  drh
	);
	EAM(std::string filename);

	// potential U(r)
	double value(const vector r);
	double value(const double* r2, std::size_t size);
	void   value(const double* r2, double* Uout, std::size_t size);
	// void   value(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// virial U'(r)/r
	matrix virial(const vector r, const double AdF);
	matrix virial(const vector* r, const double* r2, const double* AdF, std::size_t size);
	void   virial(const vector* r, const double* r2, const double* AdF, matrix* Vout, std::size_t size);
	void   virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// force 
	vector force(const vector r, const double AdF);
	vector force(const vector* r, const double* r2, const double* AdF, std::size_t size);
	void   force(const vector* r, const double* r2, const double* AdF, vector* Fout, std::size_t size);
	void   force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// charge density
	double density(const double* r2, std::size_t size);
	void   density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	// embedding energy
	void   embeddingEnergy(const double* nr, double* EF, double* dF, std::size_t size);

	//temp 
	void U_arr(double* Uout, int size);
	void dU_arr(double* Uout, int size);
	void dFemb_arr(double* Uout, int size);
	void dRho_arr(double* Uout, int size);

private:
	// functions interpolation
	void init_spline(const double* y, double* z, double H1, std::size_t N, std::size_t ib);

	std::vector<double> Femb, Rho, Vpar, dFemb, dRho, dVpar;
	int Nr1, Nrh1, N;
	double dr, drh, dr1, drh1;
};

} // potential
} // solver
} // mdcraft