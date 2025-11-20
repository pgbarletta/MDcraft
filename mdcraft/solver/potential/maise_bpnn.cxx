#include <mdcraft/solver/potential/maise_bpnn.h>

#ifdef mdcraft_ENABLE_BPNN

// #include <mdcraft/tools/timer.h>

#ifdef DEBUG
#include <mdcraft/tools/helpers.h>
using namespace mdcraft::tools;
#endif

namespace mdcraft {
namespace solver {
namespace potential {

constexpr double eV = 1.602176634;
constexpr double Na = 6.02214076;
constexpr double eV_Na = eV * Na;

MaiseBPNNPotential::MaiseBPNNPotential(
    std::string_view filename,
    Domain& domain
) : BaseTernary(0, ::mdcraft::solver::potential::TypeTernary::behler_parrinello, domain)
  , ann_(MaiseBPNN(filename))
  ,	descriptor_(MaiseBPDescriptor(filename))
{
    R1cut_ = descriptor_.Rc(); // in [nm] already!
    R2cut_ = R1cut_ * R1cut_;
}

MaiseBPNNPotential::~MaiseBPNNPotential()
{
}

void MaiseBPNNPotential::reset_natoms(std::size_t natoms) {
    tinfo_.resize(natoms);
}

void MaiseBPNNPotential::refresh(
    Atoms::iterator atomit,
    Atoms& atoms,
    const NeibsList& nlist
) {
    auto ii = std::distance(atoms.begin(), atomit);
    tinfo_[ii].g.setZero();
    tinfo_[ii].dedg.setZero();
    tinfo_[ii].dgdr.setZero();
    tinfo_[ii].dgdr_X_r.setZero();
}

void MaiseBPNNPotential::prepare(
    Atoms::iterator atomit,
    Atoms& atoms,
    const NeibsList& nlist
) {
    refresh(atomit, atoms, nlist);
    accumulate_inner_grads(atomit, atoms, nlist);
    dEdG(atomit, atoms);
}

void MaiseBPNNPotential::force(
    Atoms::iterator atomit,
    Atoms& atoms,
    const NeibsList& nlist
) {
    // Is it safe? Usually always called via stepper.make_step
    //prepare(atomit, atoms, nlist);

    auto& atom_i = *atomit;

    // we need neibs of atom here
    auto ii = std::distance(atoms.begin(), atomit);
    assert(&atoms[ii] == &atom_i);
    const auto& nlist_i = nlist[ii];

    auto in = tinfo_[ii].dgdr.dimension(0);
    assert(tinfo_[ii].dgdr.dimension(1) == 3);

    for (std::size_t j = 0; j < nlist_i.size(); j++) {
        auto jj       = nlist_i[j]; // formally ii is also considered here!

        for (std::size_t l = 0; l < in; ++l)
            for (std::size_t c = 0; c < 3; ++c)
                atom_i.f(c) += tinfo_[jj].dedg(l, 0) * tinfo_[ii].dgdr(l, c, jj);  
    }


#ifdef DEBUG
    std::cout << "force on atom ii [eV/nm] = " << ii << std::endl;
    for (std::size_t c = 0; c < 3; ++c)
        print(atom_i.f(c) / 10.0 / eV_Na);
    std::cout << "its coords:" << std::endl;
    print(atom_i.r);
    std::cout << "Energy of atom ii [eV] = " << atom_i.Ep / 10.0 / eV_Na << std::endl;
#endif
}

void MaiseBPNNPotential::virial(
    Atoms::iterator atomit,
    Atoms& atoms,
    const NeibsList& nlist
) {
    auto& atom_i = *atomit;

    // we need neibs of atom here
    auto ii = std::distance(atoms.begin(), atomit);
    assert(&atoms[ii] == &atom_i);
    const auto& nlist_i = nlist[ii];

    auto in = tinfo_[ii].dgdr.dimension(0);
    assert(tinfo_[ii].dgdr.dimension(1) == 3);

    for (std::size_t j = 0; j < nlist_i.size(); j++) {
        auto jj       = nlist_i[j]; // formally ii is also considered here!

        for (std::size_t l = 0; l < in; ++l)
            for (std::size_t c = 0; c < 3; ++c)
                atom_i.V(c, c) -= tinfo_[jj].dedg(l, 0) * tinfo_[ii].dgdr_X_r(l, c, jj);
    }
}

void MaiseBPNNPotential::accumulate_inner_grads(
    Atoms::iterator atomit,
    Atoms& atoms,
    const NeibsList& nlist
) {
    auto& atom_i = *atomit;

    // assumed atoms and neibs is the same object!
    auto ii = std::distance(atoms.begin(), atomit); // global index of atom i (this thread)
    assert(&atoms[ii] == &atom_i);
    auto& nlist_i = nlist[ii];

    #pragma ivdep
    #pragma vector aligned
    for (std::size_t j = 0, szi = nlist_i.size(); j < szi; j++) {
        auto jj       = nlist_i[j]; // global index of j-th neib of atom i
        auto& atom_j  = atoms[jj];
        auto& nlist_j = nlist[jj];

        auto kindDiff_ij = (atom_i.kind - atom_j.kind + ann_.nKinds()) % ann_.nKinds();

        vector r_ij = atom_j.r - atom_i.r;
        r_ij = domain_.shortest(r_ij);
        double const r2ij  = r_ij.squaredNorm();
        if (r2ij < 1e-15) continue;

        if (r2ij < R2cut_)
        {
            double const rij   = std::sqrt(r2ij);
            vector const dr_ij = r_ij / rij;

            // first do radial part (2 particles)
            for (std::size_t ir = 0; ir < descriptor_.numRadial(); ++ir) {
                auto pos = ir + descriptor_.shiftRadial()[kindDiff_ij];
                tinfo_[ii].g(pos, 0) += descriptor_.radial_function(ir, rij);
                for (std::size_t c = 0; c < 3; ++c) {
                    // initially tinfo_[ii].dgdr(pos, c, ii, jj) = 0 for all fixed pair of ii and jj
                    auto add = descriptor_.radial_function_d(ir, rij) * dr_ij(c);
                    tinfo_[ii].dgdr(pos, c, jj) += add;
                    tinfo_[ii].dgdr(pos, c, ii) += add;
                    tinfo_[ii].dgdr_X_r(pos, c, jj) -= add * r_ij(c);
                }
            }

            //then angular part (3 particles)
            // 1. neibs of neib j (j-th is in the center of an angle)
            #pragma ivdep
            #pragma vector aligned
            for (std::size_t k = 0, szj = nlist_j.size(); k < szj; k++) {
                auto kk      = nlist_j[k];
                auto& atom_k = atoms[kk];

                vector r_jk  = atom_k.r - atom_j.r;
                r_jk = domain_.shortest(r_jk);
                double const r2jk  = r_jk.squaredNorm();
                if (r2jk < 1e-15) continue;
                double const rjk   = std::sqrt(r2jk);
                vector const dr_jk = r_jk / rjk;

                auto kindDiff_jk = (atom_k.kind - atom_j.kind + ann_.nKinds()) % ann_.nKinds();

                if (r2jk < R2cut_) { // WRONG FOR G5 type
                    // the third triangle edge (directed from k to i atom periodic images)
                    auto r_ik = r_ij + r_jk;
#ifdef DEBUG
                    vector r_ik_ = atom_k.r - atom_i.r;
                    r_ik_ = domain_.shortest(r_ik_);
                    assert(r_ik_ == r_ik);
#endif
                    double const r2ik  = r_ik.squaredNorm();
                    if (r2ik < 1e-15) continue; // degenerate angle, atom_i == atom_k

                    if (r2ik < R2cut_) {
                        double const rik   = std::sqrt(r2ik);
                        vector const dr_ik = r_ik / rik;

                        double const cos_j = -dr_ij.dot(dr_jk);
                        auto dcos_j = [&] (unsigned short c) {
                            return r_ij(c)*cos_j/rij/rij + r_jk(c)/rij/rjk;
                        };

                        auto kindDiff_ijk = kindDiff_ij*kindDiff_ij + kindDiff_jk*kindDiff_jk;

                        /* here maximum perfomance drop! */
                        for (std::size_t ia = 0; ia < descriptor_.numAngular(); ++ia) {
                            auto pos = ia + descriptor_.shiftAngular()[kindDiff_ijk];
                            for (unsigned short c = 0; c < 3; ++c) {
                                tinfo_[ii].dgdr(pos, c, jj)
                                    -= 2.0 * descriptor_.angular_function_d(ia, rij, rik, rjk, cos_j, dcos_j(c), -dr_ij(c), -dr_ik(c));
                            }
                        }
                    }
                }

            } // angular part 1

            // 2. double run over i-th neibs (i-th is in the center of an angle)
            #pragma ivdep
            #pragma vector aligned
            for (std::size_t j1 = j, szi = nlist_i.size(); j1 < szi; j1++) {
                auto kk      = nlist_i[j1];
                auto& atom_k = atoms[kk];

                vector r_jk  = atom_k.r - atom_j.r;
                r_jk = domain_.shortest(r_jk);
                double const r2jk  = r_jk.squaredNorm();
                if (r2jk < 1e-15) continue;
                double const rjk   = std::sqrt(r2jk);
                vector const dr_jk = r_jk / rjk;

                auto kindDiff_ik = (atom_k.kind - atom_i.kind + ann_.nKinds()) % ann_.nKinds();

                if (r2jk < R2cut_) { // WRONG FOR G5 type
                    // the third triangle edge (directed from k to i atom periodic images)
                    auto r_ik = r_ij + r_jk;
#ifdef DEBUG
                    vector r_ik_ = atom_k.r - atom_i.r;
                    r_ik_ = domain_.shortest(r_ik_);
                    assert(r_ik_ == r_ik);
#endif
                    double const r2ik  = r_ik.squaredNorm();
                    if (r2ik < 1e-15) continue; // degenerate angle, atom_i == atom_k

                    if (r2ik < R2cut_) {
                        double const rik   = std::sqrt(r2ik);
                        vector const dr_ik = r_ik / rik;

                        double const cos_i =  dr_ij.dot(dr_ik);
                        auto dcos_i = [&] (unsigned short c) {
                            return cos_i * ( dr_ij(c)/rij + dr_ik(c)/rik ) - dr_ij(c)/rik - dr_ik(c)/rij;
                        };

                        auto kindDiff_ijk = kindDiff_ij*kindDiff_ij + kindDiff_ik*kindDiff_ik;

                        for (std::size_t ia = 0; ia < descriptor_.numAngular(); ++ia) {
                            auto pos = ia + descriptor_.shiftAngular()[kindDiff_ijk];
                            tinfo_[ii].g(pos, 0) += descriptor_.angular_function(ia, rij, rik, rjk, cos_i);
                            for (unsigned short c = 0; c < 3; ++c) {
                                tinfo_[ii].dgdr(pos, c, ii)
                                    -= 2.0 * descriptor_.angular_function_d(ia, rij, rik, rjk, cos_i, dcos_i(c), -dr_ij(c), -dr_ik(c));
                                tinfo_[ii].dgdr_X_r(pos, c, ii)
                                    -= 2.0 * descriptor_.angular_function_stress(ia, rij, rik, rjk, cos_i, dcos_i(c), -dr_ij(c), -dr_ik(c), r_ij(c), r_ik(c));
                            }
                        }
                    }
                }

            } // angular part 2
        }
    }
}

// TODO: just evaluate energy without gradients
void MaiseBPNNPotential::dEdG(Atoms::iterator atomit, Atoms& atoms) {
    auto& atom_i = *atomit;
    auto ii = std::distance(atoms.begin(), atomit);

    auto& g = tinfo_[ii].g;
    assert(g.size() <= INPUT_LAYER);
    //assert(std::any_of(g.data(), g.data() + g.size(), [] (auto i) { return i != 0.0; }));
    assert(std::all_of(g.data() + INPUT_LAYER, g.data() + g.size(), [] (auto i) { return i == 0.0; }));

    Eigen::DSizes<Eigen::DenseIndex, 2> offsets(0, 0);
    Eigen::DSizes<Eigen::DenseIndex, 2> extents(ann_.neurons(0), 1);
    Tens2<double> in = g.slice(offsets, extents);
    assert(in.size() == ann_.neurons(0));

    auto  kind = atom_i.kind; // different ANN in acc. w. kind

    auto& b_01 =  ann_.biases(kind, 0);
    auto& b_12 =  ann_.biases(kind, 1);
    auto& b_23 =  ann_.biases(kind, 2);

    auto& w_01 = ann_.weights(kind, 0);
    auto& w_12 = ann_.weights(kind, 1);
    auto& w_23 = ann_.weights(kind, 2);

    // explicitly write code for 4 layers
    Tens2<ADiffd> dB_01 = mdcraft::bpnn::autodiff::prepare( b_01 );
    Tens2<ADiffd> dW_01 = mdcraft::bpnn::autodiff::prepare( w_01 );
    Tens2<ADiffd> dX_01 = mdcraft::bpnn::autodiff::prepare(
                   Eigen::TensorMap<Tens2<double>>{in.data(), in.size(), 1l}, 0, in.size());
    Tens2<ADiffd> dY_01 = Layer<ADiffd>(dX_01, dB_01, dW_01, tanh_af<ADiffd>);

    Tens2<ADiffd> dB_12 = mdcraft::bpnn::autodiff::prepare( b_12 );
    Tens2<ADiffd> dW_12 = mdcraft::bpnn::autodiff::prepare( w_12 );
    Tens2<ADiffd> dY_12 = Layer<ADiffd>(dY_01, dB_12, dW_12, tanh_af<ADiffd>);

    Tens2<ADiffd> dB_23 = mdcraft::bpnn::autodiff::prepare( b_23 );
    Tens2<ADiffd> dW_23 = mdcraft::bpnn::autodiff::prepare( w_23 );
    Tens2<ADiffd> dY_23 = Layer<ADiffd>(dY_12, dB_23, dW_23, linear_af<ADiffd>); // linear activation here

    assert(dY_23.dimension(0) == dY_23.dimension(1) == 1);
    auto E = static_cast<Eigen::Tensor<ADiffd, 0>>(dY_23.sum())(0);

    atom_i.Ep = E.value() * 10.0 * eV_Na; // from [eV] to [kJ/mol];
    Tens2<double> out = grads<double>(E, dX_01) * 10.0 * eV_Na;
    Eigen::array<Eigen::IndexPair<long int>, 2> padding;
    padding[0] = {0, INPUT_LAYER - ann_.neurons(0)};
    padding[1] = {0, 0};
    tinfo_[ii].dedg = out.pad(padding);

    return;
}

} // potential
} // solver
} // mdcraft

#endif