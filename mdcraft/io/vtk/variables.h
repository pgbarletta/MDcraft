#pragma once

#include <functional>

#include <mdcraft/data/atom.h>
#include <mdcraft/io/vtk/vtk_type.h>

namespace mdcraft::io {

/// @brief \~russian Тип функции для записи переменных, позволяет инициализировать
/// функцию записи переменных через лямбда функцию. Далее пример использования.
/// Позволяет сократить размер выходного файла.
///  \~english Variable writing function type, enables initialization
/// of the variable writing function via a lambda expression. Usage example below.
/// Helps reduce the size of the output file.
/// @code
/// WriteFunction<float> ev_energy =
///     [](Atom& atom, float* out) {
///         out[0] = static_cast<float>(atom.Ep / 1.6e-19);
///     };
/// @endcode
template <typename T>
using WriteAtom = std::function<void(data::Atom&, T*)>;

/** \brief
    \~russian Класс для записи переменных в VTU файл, каждой переменной для
              записи должен соответствовать экземпляр Variable.
    \~english Class for writing variables to VTU file, each variable
              to be written must correspond to a Variable instance.
    \~
*/
class Variable {
public:
    /// @brief Creating a variable without a description is forbidden
    Variable() = delete;

    /// @brief Get a variable from a preinstalled set
    static Variable Default(const std::string &name);

    /// @brief Создать переменную из функции вида Atom& -> T.
    /// @code
    /// Variable var("tag", [](Atom& atom) -> unsigned int { return atom.uid; });
    /// @endcode
    template<class Func, typename = std::enable_if_t<std::is_invocable_v<Func, data::Atom&>>>
    Variable(const std::string& name, Func&& func);

    /// @brief Create a variable with a full description
    ///
    /// The following code adds a vector variable "momentum" that allows
    /// recording two momentum components in Float32 format.
    /// @code
    /// Variable fd("momentum", 2,
    ///     WriteFunction<float>([](Atom& atom, float* out) {
    ///         out[0] = static_cast<float>(atom.mass * atom.velocity.x);
    ///         out[1] = static_cast<float>(atom.mass * atom.velocity.y);
    ///     }));
    /// @endcode
    template<class T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    Variable(const std::string& name, int n_components, const WriteAtom<T> &func);

    /// @brief Variable name
    std::string name() const { return m_name; }

    /// @brief Variable type
    VtkType type() const { return m_type; }

    /// @brief Number of components for a vector variable
    int n_components() const { return m_n_components; }

    /// @brief If a variable is a scalar
    bool is_scalar() const { return m_n_components < 2; }

    /// @brief Analog of sizeof for a variable
    size_t size() const { return m_n_components * m_type.size(); }

    /// @brief Write a variable into a buffer
    void write(data::Atom &atom, std::byte *out) const {
        m_write(atom, out);
    }

private:
    std::string m_name;  ///< Variable name
    VtkType m_type;      ///< Variable type
    int m_n_components;  ///< Components number (for vector/matrix)

    WriteAtom<std::byte> m_write = nullptr; ///< Writing function
};

/** \brief
    \~russian Список переменных для записи в VTU файл
    \~english A list of variables to be written into a VTU file
    \~
*/
class Variables {
public:
    /// @brief An empty list of variables
    Variables() = default;

    /// @brief Create a list from one name like Variables = "r"
    explicit Variables(const std::string &name);

    /// @brief Create a list from a set of variables' names
    /// Example. Variables list = {"r", "f"}
    Variables(const std::vector<std::string> &names);

    /// @brief Create a list from a set of variables' names
    /// Пример. Variables list = {"r", "f"}
    Variables(std::initializer_list<const char*> names);

    /// @brief Add a variable with a name
    void operator+=(const std::string &name);

    /// @brief Add a set of variables from their names into the list
    void operator+=(const std::vector<std::string> &names);

    /// @brief Add a set of variables from their names into the list
    void operator+=(std::initializer_list<const char*> names);

    /// @brief Add a variable from its Variable ctor parameters
    template <typename... Args>
    void append(Args&&... args) {
        m_list.emplace_back(std::forward<Args>(args)...);
    }

    /// @brief Clear the variables' list
    void clear() { m_list.clear(); }

    /// @brief The number of variables
    size_t size() const { return m_list.size(); }

    /// @brief An access to the list of variables
    const std::vector<Variable> &list() const { return m_list; }

private:
    std::vector<Variable> m_list;  ///< A list of variables
};

template<class Func, typename>
Variable::Variable(const std::string& name, Func&& func) {
    using T = std::result_of_t<Func(data::Atom&)>;

    static_assert(std::is_same_v<T, data::vector> || std::is_same_v<T, data::matrix> ||
        VtkType::get<T>() != VtkType::Undefined, "Only arithmetic type, data::vector or data::matrix");

    if constexpr (VtkType::get<T>() != VtkType::Undefined) {
        m_name = name;
        m_n_components = 1;
        m_type = VtkType::get<T>();
        m_write = [func=std::forward<Func>(func)](data::Atom& atom, std::byte* arg) {
            auto out = reinterpret_cast<T*>(arg);
            out[0] = func(atom);
        };
        return;
    }
    if constexpr (std::is_same_v<T, data::vector>) {
        m_name = name;
        m_n_components = 3;
        m_type = VtkType::Float64;
        m_write = [func=std::forward<Func>(func)](data::Atom& atom, std::byte *arg) {
            auto out = reinterpret_cast<data::vector*>(arg);
            out[0] = func(atom);
        };
        return;
    }
    if constexpr (std::is_same_v<T, data::matrix>) {
        m_name = name;
        m_n_components = 6;
        m_type = VtkType::Float64;
        m_write = [func=std::forward<Func>(func)](data::Atom& atom, std::byte *arg) {
            data::matrix M = func(atom);
            auto out = reinterpret_cast<double*>(arg);
            out[0] = M(0, 0);
            out[1] = M(1, 1);
            out[2] = M(2, 2);
            out[3] = M(0, 1);
            out[4] = M(1, 2);
            out[5] = M(0, 2);
        };
        return;
    }
}

template<class T, typename>
Variable::Variable(const std::string& name, int n_components, const WriteAtom<T> &func) {
    static_assert(VtkType::get<T>() != VtkType::Undefined);

    m_name = name;
    m_type = VtkType::get<T>();
    m_n_components = n_components;
    m_write = [func](data::Atom &atom, std::byte *out) {
        func(atom, reinterpret_cast<T*>(out));
    };
}

} // namespace mdcraft::io