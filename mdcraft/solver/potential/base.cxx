#include <mdcraft/solver/potential/base.h>

namespace mdcraft { 
namespace solver { 
namespace potential {


Base::Base(double Rcutoff, 
	::mdcraft::solver::potential::type type
) : type(type) {
	R1cut = Rcutoff;
	R2cut = R1cut * R1cut;
}

double Base::rcut() { return R1cut; }

double Base::value(const vector r) { return 0.0; }
double Base::value(const double* r2, std::size_t size) { return 0.0; }
void   Base::value(const double* r2, double* Uout, std::size_t size) { return; }
void   Base::value(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) { return; }

matrix Base::virial(const vector r) { return matrix::Zero(); };
matrix Base::virial(const vector* r, const double* r2, std::size_t size) { return matrix::Zero(); }
void   Base::virial(const vector* r, const double* r2, matrix* Vout, std::size_t size) { return; }
void   Base::virial(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) { return; }

vector Base::force(const vector r) { return vector::Zero(); };
vector Base::force(const vector* r, const double* r2, std::size_t size) { return vector::Zero(); }
vector Base::force(const vector* r, const double* r2, const double* AdF, std::size_t size) { return vector::Zero(); }
void   Base::force(const vector* r, vector* Fout, std::size_t size) { return; }
void   Base::force(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) { return; }

double Base::density(const double* r2, std::size_t size) { return 0.0; }
void   Base::density(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) { return; }

void   Base::embeddingEnergy(const double* nr, double* EF, double* dF, std::size_t size){ return; }

Base::~Base() {}

BaseManyBody::BaseManyBody(
	double Rcutoff,
	::mdcraft::solver::potential::TypeManyBody type,
	Domain& domain
) : R1cut_(Rcutoff)
  , R2cut_(Rcutoff*Rcutoff)
  , type_(type)
  , domain_(domain)
{}

BaseTernary::BaseTernary(
	double Rcutoff,
	::mdcraft::solver::potential::TypeTernary type,
	Domain& domain
) : R1cut_(Rcutoff)
  , R2cut_(Rcutoff*Rcutoff)
  , type_(type)
  , domain_(domain)
{}

BaseManyBody::~BaseManyBody()
{
}

double 				 BaseManyBody::rcut()   const { return R1cut_; }
unsigned short BaseManyBody::mytype() const { return type_;  }

// void BaseTernary::prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) { return; }
void BaseManyBody::force(   Atoms& atoms, Atoms& neibs, NeibsList& nlist) { return; }
// void BaseTernary::virial(  Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) { return; }

BaseTernary::~BaseTernary()
{
}

double 				 BaseTernary::rcut()   const { return R1cut_; }
unsigned short BaseTernary::mytype() const { return type_;  }

void BaseTernary::prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) { return; }
void BaseTernary::force(   Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) { return; }
void BaseTernary::virial(  Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) { return; }

void BaseTernary::reset_natoms(std::size_t natoms) { return; }

} // potential
} // solver
} // mdcraft
