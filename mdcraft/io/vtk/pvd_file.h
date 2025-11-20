#pragma once

#include <mdcraft/configuration.h>
#include <mdcraft/tools/network.h>
#include <mdcraft/io/vtk/vtu_file.h>

namespace mdcraft::io {

/** \brief
    \~russian Запись серии VTU-файлов и соответствующего PVD-файла
    \~english Recording a series of VTU-files and the corresponding PVD-file
    \~
*/
class PvdFile {
public:
    Variables variables;   ///< A list of variables to write

    /// @brief Opens a PVD-file, writes header
    /// @param filename A short name without an extension
    /// @param directory A relative or absolute path
    explicit PvdFile(const std::string& filename,
                     const std::string& directory = "output");

#ifdef mdcraft_ENABLE_MPI
    /// @brief Opens a PVD-file, writes header
    /// @param network The net used in parallel writing of several files
    /// @param filename A short name without an extension
    /// @param directory A relative or absolute path
    PvdFile(tools::Network network, const std::string& filename,
            const std::string& directory = "output");
#endif

    /// @brief Writing to a storage (or part, in distributed case)
    /// into a VTU-file (a set of VTU), then renew the PVD-file.
    /// VtuFile::save function is used
    void save(data::Atoms& atoms, double timestep);

private:
    PvdFile();

    void open(const std::string& filename, const std::string& _directory, bool distributed);

    std::string get_filename() const;

    void update_pvd(double timestep);

    bool           m_open;         ///< If the PVD-file is opened
    bool           m_distributed;  ///< If a PVD file is shared in MPI
    std::string    m_filename;     ///< A name without an extension
    std::string    m_fullname;     ///< An absolute name without an extension
    std::streamoff m_pos;          ///< A pointer to write in file
    std::size_t    m_counter;      ///< A counter for timesteps being written

    tools::Network m_net;
};

} // namespace mdcraft::io