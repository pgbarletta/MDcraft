#include <mdcraft/configuration.h>
#include <mdcraft/io/vtk/variables.h>
#ifdef mdcraft_ENABLE_MPI
#include <mpi.h>
#endif

namespace mdcraft::io {

using data::Atom;
using data::vector;

class Variable Variable::Default(const std::string &name) {
    if (name == "r") {
        return {name, [](Atom& atom) -> vector { return atom.r; }};
    }
    if (name == "v") {
        return {name, [](Atom& atom) -> vector { return atom.v; }};
    }
    if (name == "f") {
        return {name, [](Atom& atom) -> vector { return atom.f; }};
    }
    if (name == "m") {
        return {name, [](Atom& atom) -> double { return atom.m; }};
    }
    if (name == "rcut") {
        return {name, [](Atom& atom) -> double { return atom.rcut; }};
    }
    if (name == "rns") {
        return {name, [](Atom& atom) -> double { return atom.rns; }};
    }
    if (name == "Ep") {
        return {name, [](Atom& atom) -> double { return atom.Ep; }};
    }
    if (name == "Em") {
        return {name, [](Atom& atom) -> double { return atom.Em; }};
    }
    if (name == "nc") {
        return {name, [](Atom& atom) -> double { return atom.nc; }};
    }
    if (name == "bc") {
        return {name, [](Atom& atom) -> double { return atom.bc; }};
    }
    if (name == "kind") {
        return {name, [](Atom& atom) -> unsigned short { return atom.kind; }};
    }
    if (name == "uid") {
        return {name, [](Atom& atom) -> unsigned short { return atom.uid; }};
    }
    if (name == "tag") {
        return {name, [](Atom& atom) -> unsigned int { return atom.tag; }};
    }
    if (name == "rank") {
#ifndef mdcraft_ENABLE_MPI
        return {name, [](Atom& atom) -> int { return 0; }};
#else
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return {name, [rank](Atom& ) -> int { return rank; }};
#endif
    }

    throw std::runtime_error("Unknown default variable '" + std::string(name) + "'");
}

Variables::Variables(const std::string& name) {
    operator+=(name);
}

Variables::Variables(const std::vector<std::string>& names) {
    for (const std::string& name: names) { operator+=(name); }
}

Variables::Variables(std::initializer_list<const char*> names) {
    for (const char *name: names) { operator+=(name); }
}

void Variables::operator+=(const std::string& name) {
    m_list.emplace_back(Variable::Default(name));
}

void Variables::operator+=(const std::vector<std::string> &names) {
    for (auto& name: names) { operator+=(name); }
}

void Variables::operator+=(std::initializer_list<const char *> names) {
    for (auto& name: names) { operator+=(name); }
}

} // namespace mdcraft::io