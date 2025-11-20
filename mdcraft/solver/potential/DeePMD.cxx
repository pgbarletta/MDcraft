#include <mdcraft/solver/potential/DeePMD.h>

#ifdef mdcraft_ENABLE_DeePMD

namespace mdcraft {
namespace solver {
namespace potential {

constexpr double eV = 1.602176634;
constexpr double Na = 6.02214076;
constexpr double eV_Na = eV * Na;


// taken explicitly from DeePMD
std::string get_file_content(const std::string &model) {
  int myrank = 0, root = 0;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int nchar = 0;
  std::string file_content;
  if (myrank == root) {
    deepmd_compat::read_file_to_string(model, file_content);
    nchar = file_content.size();
  }
  // MPI_Bcast(&nchar, 1, MPI_INT, root, MPI_COMM_WORLD);
  char *buff = (char *)malloc(sizeof(char) * nchar);
  if (myrank == root) {
    memcpy(buff, file_content.c_str(), sizeof(char) * nchar);
  }
  // MPI_Bcast(buff, nchar, MPI_CHAR, root, MPI_COMM_WORLD);
  file_content.resize(nchar);
  for (unsigned ii = 0; ii < nchar; ++ii) {
    file_content[ii] = buff[ii];
  }
  free(buff);
  return file_content;
}


DeepModelingPotential::DeepModelingPotential(
    const std::string& filename,
    Domain& domain
) : BaseManyBody(0, ::mdcraft::solver::potential::TypeManyBody::deep_modeling, domain)
{
    std::vector<std::string> models;
    numb_models = 1;//models.size();
    if (numb_models == 1) {
        try {
            deep_pot.init(filename, 0, get_file_content(filename));
        } catch (deepmd_compat::deepmd_exception &e) {
            std::cout << "[DeepModelingPotential] ERROR : " << e.what() << std::endl;
        }
        cutoff = deep_pot.cutoff() * 0.1; // [A] -> [nm]
        numb_types = deep_pot.numb_types();
        numb_types_spin = deep_pot.numb_types_spin();
        dim_fparam = deep_pot.dim_fparam();
        dim_aparam = deep_pot.dim_aparam();

        fparam.clear();
        aparam.clear();

        // std::cout << "cutoff = " << cutoff << '\n';
        // std::cout << "numb_types = " << numb_types << '\n';
        // std::cout << "numb_types_spin = " << numb_types_spin << '\n';
        // std::cout << "dim_fparam = " << dim_fparam << '\n';
        // std::cout << "dim_aparam = " << dim_aparam << '\n';
    }

    R1cut_ = cutoff; // in [nm] already!
    R2cut_ = R1cut_ * R1cut_;
}

void make_uniform_aparam(std::vector<double> &daparam,
                         const std::vector<double> &aparam,
                         const int &nlocal) {
  unsigned dim_aparam = aparam.size();
  daparam.resize(static_cast<size_t>(dim_aparam) * nlocal);
  for (int ii = 0; ii < nlocal; ++ii) {
    for (int jj = 0; jj < dim_aparam; ++jj) {
      daparam[ii * dim_aparam + jj] = aparam[jj];
    }
  }
}


DeepModelingPotential::~DeepModelingPotential()
{
}

void DeepModelingPotential::force(Atoms& atoms, Atoms& neibs, NeibsList& nlist)
{
    auto nall = atoms.size();
    auto nlocal = nall;

    std::vector<int> dtype(nall);
    for (int ii = 0; ii < nall; ++ii) {
        dtype[ii] = atoms[ii].kind;
    }

    double dener(0);
    std::vector<double> dforce(nall * 3);
    std::vector<double> dvirial(9, 0);
    std::vector<double> dcoord(nall * 3, 0.);
    std::vector<double> dbox(9, 0);
    std::vector<double> daparam;

    std::vector<double> dener_per_atom(nall, 0.);
    std::vector<double> dvirial_per_atom(nall * 9, 0.);

    // get box
    dbox[0] = domain_.xmax() * 10.;  // xx
    dbox[4] = domain_.ymax() * 10.;  // yy
    dbox[8] = domain_.zmax() * 10.;  // zz

    // get coord
    for (int ii = 0; ii < nall; ++ii) {
        auto& atom_i = atoms[ii];
        dcoord[ii * 3 + 0] = (atom_i.r(0) - domain_.xmin()) * 10.; 
        dcoord[ii * 3 + 1] = (atom_i.r(1) - domain_.ymin()) * 10.; 
        dcoord[ii * 3 + 2] = (atom_i.r(2) - domain_.zmin()) * 10.; 
    }

    if (aparam.size() > 0) {
        // uniform aparam
        make_uniform_aparam(daparam, aparam, nlocal);
    }

    /*  Leave it here: via NeibList

    int nghost = 0;
    int ago = 0;

    auto inum = atoms.size();
    std::vector<int> ilist(inum);
    for (size_t i = 0; i < inum; ++i) ilist[i] = i;

    std::vector<int> numneigh(inum);
    std::transform(nlist.begin(), nlist.end(), numneigh.begin(), [] (const auto& l) { return l.size() - 1; });

    int** firstneigh = new int*[inum];

    for (int i = 0; i < inum; ++i) {
        int jj = 0;
        int* neib_inds = new int[numneigh[i]];
        for (int j = 0; j < numneigh[i] + 1; ++j) {
            if (&atoms[nlist[i][j]] == &atoms[i]) continue;
            neib_inds[jj++] = nlist[i][j];
        }
        firstneigh[i] = neib_inds;
    }
    deepmd_compat::InputNlist list_(inum, ilist.data(), numneigh.data(), firstneigh);
    */

    try {
      deep_pot.compute(dener, dforce, dvirial, dener_per_atom, dvirial_per_atom,
                       dcoord, dtype, dbox, fparam, daparam);
      /* via NeibList */
      // deep_pot.compute(dener, dforce, dvirial, dener_per_atom, dvirial_per_atom,
                       // dcoord, dtype, dbox, nghost, list_, ago, fparam, daparam);
    } catch (deepmd_compat::deepmd_exception &e) {
      std::cout << "[DeePMD::force] ERROR : " << e.what() << std::endl;
    }

    // get force
    for (int ii = 0; ii < nall; ++ii) {
        auto& atom_i = atoms[ii];

        atom_i.Ep = dener_per_atom[ii] * 10. * eV_Na;

        atom_i.f(0) += dforce[3 * ii + 0] * 100. * eV_Na;
        atom_i.f(1) += dforce[3 * ii + 1] * 100. * eV_Na;
        atom_i.f(2) += dforce[3 * ii + 2] * 100. * eV_Na;

        atom_i.V(0,0) -= dvirial_per_atom[9 * ii + 0] * 10. * eV_Na; // xx
        atom_i.V(1,1) -= dvirial_per_atom[9 * ii + 4] * 10. * eV_Na; // yy
        atom_i.V(2,2) -= dvirial_per_atom[9 * ii + 8] * 10. * eV_Na; // zz

        atom_i.V(0,1) -= dvirial_per_atom[9 * ii + 3] * 10. * eV_Na; // xy
        atom_i.V(0,2) -= dvirial_per_atom[9 * ii + 4] * 10. * eV_Na; // xz
        atom_i.V(1,2) -= dvirial_per_atom[9 * ii + 5] * 10. * eV_Na; // yz

        atom_i.V(1,0) -= atom_i.V(0,1);
        atom_i.V(2,0) -= atom_i.V(0,2);
        atom_i.V(2,1) -= atom_i.V(1,2);
    }
}

} // potential
} // solver
} // mdcraft

#endif