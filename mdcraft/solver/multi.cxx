#include <cmath>
#include <mdcraft/solver/multi.h>

namespace mdcraft {
namespace solver {

using ::mdcraft::tools::Threads;
using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

Multi::Multi(
	Boundary&   boundary,
	Thermostat& thermostat,
	std::size_t nkinds,
	Threads&    threads
) : nkinds(nkinds),
	Pair(
		boundary,
	    thermostat,
	    threads 
    )
{
	potentials.resize(nkinds, std::vector<Potential*>(nkinds));
}

void Multi::add_potential(std::size_t i, std::size_t j, Potential& pot) {
	potentials[i][j] = &pot;
	potentials[j][i] = &pot;

	if (stages.empty()) {
		if (pot.type == potential::type::pair) {
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				forces(atoms, neibs, nlist);
			});
		}
		else if (pot.type == potential::type::eam) {
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				charge(atoms, neibs, nlist);
			});
			stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
				forces(atoms, neibs, nlist);
			});
		}
	} else if(pot.type == potential::type::eam and stages.size() == 1) {
		stages.emplace(stages.begin(), [this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
			charge(atoms, neibs, nlist);
		});
	}
}

void Multi::prepare_one(
	Atoms::iterator atomit, 
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	atom.f  = vector::Zero();
	atom.V  = matrix::Zero();
	atom.Ep = 0.0;
	atom.Em = 0.0;
	atom.nc = 0.0;
}

void Multi::charge_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
){
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();
	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist.reserve(nlist_size);
	}

	auto const atomKind = atom.kind;

	for (auto i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		auto neib = neibs[j];
		auto const neibKind = neib.kind;
		
		// if (potentials[atomKind][neibKind]->type == 
		//     potential::type::pair) continue;

		multi_nlist[neibKind].push_back(j);
	}
	
	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		potentials[atomKind][neibKind]->density(atomit, neibs, multi_nlist[neibKind]);
	}
}

void Multi::force_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();
	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist.reserve(nlist_size);
	}

	auto const atomKind = atom.kind;

	for (auto i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		auto neib = neibs[j];
		auto const neibKind = neib.kind;
		multi_nlist[neibKind].push_back(j);
	}

	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		potentials[atomKind][neibKind]->force(atomit, neibs, multi_nlist[neibKind]);
	}
}

void Multi::virial_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	auto const nlist_size = nlist.size();
	std::vector<NeibsListOne> multi_nlist(nkinds);

	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		multi_nlist.reserve(nlist_size);
	}

	auto const atomKind = atom.kind;

	for (auto i = 0; i < nlist_size; ++i) {
		auto& j = nlist[i];
		auto neib = neibs[j];
		auto const neibKind = neib.kind;
		multi_nlist[neibKind].push_back(j);
	}

	for (auto neibKind = 0; neibKind < nkinds; ++neibKind) {
		potentials[atomKind][neibKind]->virial(atomit, neibs, multi_nlist[neibKind]);
	}
}

} // solver
} // mdcraft