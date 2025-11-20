#pragma once

#include <string>
#include <iostream>

namespace mdcraft::io {

/** \brief
    \~russian Тип скаляра при сохранении в VTU, работает как enum class.
    \~english The type of scalar when saving to VTU, works as a enum class
    \~
*/
class VtkType {
public:
    enum ValueType : int {
        Undefined = 0,
        Int8  = 1, Int16  = 2, Int32  = 3, Int64  = 4,
        UInt8 = 5, UInt16 = 6, UInt32 = 7, UInt64 = 8,
        Float32 = 9, Float64 = 10
    };

    /// @brief Ctor of VtkType type = VtkType::Int32
    constexpr VtkType(ValueType value = Undefined)
        : m_value(value <= Float64 ? value : Undefined) {
    }

    /// @brief Ctor of VtkType type = "Int32"
    VtkType(const char* value);

    /// @brief Ctor of VtkType type = "Int32"
    VtkType(const std::string& value);

    /// @brief Conversion to string
    const char* to_str() const;

    /// @brief Conversion to C-style string
    operator const char*() const { return to_str(); }

    /// @brief Conversion to C++-style string
    operator std::string() const { return to_str(); }

    constexpr bool operator==(VtkType value) const {
        return m_value == value.m_value;
    }

    constexpr bool operator!=(VtkType value) const {
        return m_value != value.m_value;
    }

    /// @brief Size of the type in bytes
    size_t size() const;

    /// @brief Get VtkType from a C++ type
    /// VtkType type = VtkType::get<int>();
    template <class T>
    static constexpr VtkType get();

protected:
    ValueType m_value;
};

template <class T>
constexpr VtkType VtkType::get() {
    if constexpr (!std::is_arithmetic_v<T>) {
        return VtkType::Undefined;
    }
    if constexpr (std::is_floating_point_v<T>) {
        // floating type
        switch (sizeof(T)) {
            case 4:  return VtkType::Float32;
            case 8:  return VtkType::Float64;
            default: return VtkType::Undefined;
        }
    }
    if constexpr (std::is_signed_v<T>) {
        // signed integral types
        switch (sizeof(T)) {
            case 1:  return VtkType::Int8;
            case 2:  return VtkType::Int16;
            case 4:  return VtkType::Int32;
            case 8:  return VtkType::Int64;
            default: return VtkType::Undefined;
        }
    }
    // unsigned integral types
    switch (sizeof(T)) {
        case 1:  return VtkType::UInt8;
        case 2:  return VtkType::UInt16;
        case 4:  return VtkType::UInt32;
        case 8:  return VtkType::UInt64;
        default: return VtkType::Undefined;
    }
}

} // namespace mdcraft::io

