#include <mdcraft/solver/potential/mlip4.h>

#ifdef mdcraft_ENABLE_MLIP4

namespace mdcraft {
namespace solver {
namespace potential {

constexpr double eV = 1.602176634;
constexpr double Na = 6.02214076;
constexpr double eV_Na = eV * Na;

MLIP4Potential::MLIP4Potential(
        const std::string& pot_filename,
        Domain& domain
    ) : BaseManyBody(0, ::mdcraft::solver::potential::TypeManyBody::mlip_4, domain)
{
    try {
        std::ifstream ifs(pot_filename, std::ios::binary);
        json_io::Reader r(ifs);
        if ((r >> mtp_pot).is_error()) throw std::runtime_error("");
        R1cut_ = mtp_pot->Cutoff() * 0.1; // [A] -> [nm]
        R2cut_ = R1cut_ * R1cut_;
    }
    catch (...) {
        throw std::runtime_error("Can not read a potential from \"" + pot_filename + "\"\n");
    }
}

void MLIP4Potential::create_nbhoods(Atoms& atoms, Atoms& neibs, NeibsList& nlist)
{
	auto nat = atoms.size();
	frame_atoms.resize(nat);

	for (auto it = atoms.begin(); it != atoms.end(); ++it) {
		auto i = it - atoms.begin();
		auto& atom_i = *it;
		auto& list_i = nlist[i];

		NBhood nbh;
		nbh.my_ind_ = i;
		nbh.my_kind_ = atom_i.kind;
		nbh.my_pos = atom_i.r;

		nbh.size_ = list_i.size() - 1;

		nbh.inds = std::vector<int>(nbh.size_);
		nbh.kinds = std::vector<int>(nbh.size_);
		nbh.rel_pos_ = std::vector<vector>(nbh.size_, vector{});

		int jj = 0;
		for (auto j : list_i) {
			if (&neibs[j] == &atom_i) continue;

			nbh.inds[jj] = j; // formally: neibs.begin() + j - atoms.begin(); //
			nbh.kinds[jj] = neibs[j].kind;

            auto dr = neibs[j].r - atom_i.r;
			nbh.rel_pos_[jj] = domain_.shortest(dr) * 10.; // from [nm] to [A]

			jj++;
		}

		frame_atoms[i] = nbh;
	}
}

void MLIP4Potential::force(Atoms& atoms, Atoms& neibs, NeibsList& nlist)
{
    auto nat = atoms.size();

    // generate frame_atoms from atoms
    create_nbhoods(atoms, neibs, nlist);

    std::vector<vector> forces(nat, vector{});

    for (auto it = frame_atoms.begin(); it != frame_atoms.end(); ++it)
    {
    	const auto& nbh = *it;

        double energy{0.0};

    	static_cast<AnyPairDescriptorPot*>(mtp_pot.get())
    		->AccumulateSiteEnergyDers_template< NBhood, Eigen::Matrix<double, 3, 1>*  >(nbh, energy, forces.data());

    	atoms[nbh.my_ind()].Ep += energy * 10. * eV_Na; // ???
    }

    // after that we have all the energies and all the forces calculated

    // extract forces
    for (auto it = frame_atoms.begin(); it != frame_atoms.end(); ++it)
    {
        const auto& nbh = *it;
        auto& atom_i = atoms[nbh.my_ind()];

        atom_i.f = forces[nbh.my_ind()] * 100. * eV_Na;
    }

    // extract virials
    for (auto it = frame_atoms.begin(); it != frame_atoms.end(); ++it)
    {
        const auto& nbh = *it;
        auto& atom_i = atoms[nbh.my_ind()];

        // loop over my nbhood
        for (auto j = 0; j < nbh.size(); ++j) {
            auto& atom_j = neibs[nbh.ind(j)];

            double fx = atom_j.f(0);
            double fy = atom_j.f(1);
            double fz = atom_j.f(2);

            double dx = nbh.rel_pos(j)[0] * 0.1;
            double dy = nbh.rel_pos(j)[1] * 0.1;
            double dz = nbh.rel_pos(j)[2] * 0.1;

            atom_i.V(0,0) -= fx * dx; // xx
            atom_i.V(1,1) -= fy * dy; // yy
            atom_i.V(2,2) -= fz * dz; // zz

            atom_i.V(0,1) -= 0.5 * (fy * dx + fx * dy); // xy
            atom_i.V(0,2) -= 0.5 * (fz * dx + fx * dz); // xz
            atom_i.V(1,2) -= 0.5 * (fz * dy + fy * dz); // yz

            atom_i.V(1,0) -= atom_i.V(0,1);
            atom_i.V(2,0) -= atom_i.V(0,2);
            atom_i.V(2,1) -= atom_i.V(1,2);
        }
    }
}

MLIP4Potential::~MLIP4Potential() {}

} // potential
} // solver
} // mdcraft

#endif