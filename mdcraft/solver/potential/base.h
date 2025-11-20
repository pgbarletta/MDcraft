#pragma once

#include <mdcraft/data/vector.h>
#include <mdcraft/data/matrix.h>
#include <mdcraft/data/atom.h>

#include <mdcraft/lattice/domain.h>
#include <mdcraft/neibs/verlet-list.h>


namespace mdcraft { 
namespace solver { 
namespace potential {

enum type : unsigned short {
	undefined,
	pair,
	eam,
	mtpr,
	ternary
};

using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;
using ::mdcraft::data::Atoms;

using NeibsList    = ::mdcraft::neibs::List;
using NeibsListOne = ::mdcraft::neibs::ListOne;

class Base {
public:
	Base(double Rcutoff, 
		::mdcraft::solver::potential::type type = 
		::mdcraft::solver::potential::type::undefined
	);
	// cutoff radius
	double rcut();
	// potential U(r)
	virtual double value(const vector r);// { return 0.0; }
	virtual double value(const double* r2, std::size_t size);// { return 0.0; }
	virtual void   value(const double* r2, double* Uout, std::size_t size);// { return; }
	virtual void   value(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);// { return; }
	// virial U'(r)/r
	virtual matrix virial(const vector r);// { return matrix::Zero(); };
	virtual matrix virial(const vector* r, const double* r2, std::size_t size);// { return matrix::Zero(); }
	virtual void   virial(const vector* r, const double* r2, matrix* Vout, std::size_t size);// { return; }
	virtual void   virial(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);// { return; }
	// force
	virtual vector force(const vector r);// { return vector::Zero(); };
	virtual vector force(const vector* r, const double* r2, std::size_t size);// { return vector::Zero(); }
	virtual vector force(const vector* r, const double* r2, const double* AdF, std::size_t size);// { return vector::Zero(); }
	virtual void   force(const vector* r, vector* Fout, std::size_t size);// { return; }
	virtual void   force(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);// { return; }

	virtual double density(const double* r2, std::size_t size);// { return 0.0; }
	virtual void   density(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist);// { return; }

	virtual void   embeddingEnergy(const double* nr, double* EF, double* dF, std::size_t size);// { return; }

	virtual ~Base();// {}

	const unsigned short type;
protected:
	double R1cut, R2cut;
};

enum TypeManyBody : unsigned short {
	default_many_body,
	deep_modeling,
	mlip_4
};

enum TypeTernary : unsigned short {
	default_ternary,
	behler_parrinello
};

using ::mdcraft::lattice::Domain;
using ::mdcraft::lattice::dummy_domain;

class BaseManyBody {
public:
	BaseManyBody(
		double Rcutoff,
		::mdcraft::solver::potential::TypeManyBody type = 
		::mdcraft::solver::potential::TypeManyBody::default_many_body,
                Domain& domain = dummy_domain
        );

        double                            rcut()   const;
        unsigned short            mytype() const;

	//virtual void          prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist);
        virtual void          force(   Atoms& atoms, Atoms& neibs, NeibsList& nlist);
        // virtual void          virial(  Atoms& atoms, Atoms& neibs, const NeibsList& nlist);

        virtual ~BaseManyBody();

protected:
        double R1cut_, R2cut_;
        unsigned short type_;
        Domain domain_;
};

class BaseTernary {
public:
	BaseTernary(
		double Rcutoff, 
		::mdcraft::solver::potential::TypeTernary type = 
		::mdcraft::solver::potential::TypeTernary::default_ternary,
		Domain& domain = dummy_domain
	);

	double 				  rcut()   const;
	unsigned short 		  mytype() const;

	virtual void          prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist);
	virtual void          force(   Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist);
	virtual void          virial(  Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist);

	/* a helper to be used in MaiseeBPNNPotential
	   and maybe others that depend on the size of Atoms they work */
	virtual void 		  reset_natoms(std::size_t natoms);

	virtual ~BaseTernary();

protected:
	double R1cut_, R2cut_;
	unsigned short type_;
	Domain domain_;
};

} // potential
} // solver
} // mdcraft
