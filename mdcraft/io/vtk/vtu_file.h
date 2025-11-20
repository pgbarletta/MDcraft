#pragma once

#include <mdcraft/io/vtk/variables.h>

namespace mdcraft::io {

/** \brief
    \~russian Запись частиц в VTU-файл
    \~english Writing particles into a VTU-file
    \~
*/
class VtuFile {
public:
    std::string filename;   ///< Full name of the file

    /// @brief Create file
    explicit VtuFile(std::string filename);

    /// @brief Write to file 
    void save(data::Atoms &atoms, const Variables& variables = {}) const;

    /// @brief Static method to write to file
    static void save(std::string filename, data::Atoms &atoms,
                     const Variables &variables = {});
};

} // namespace mdcraft::io