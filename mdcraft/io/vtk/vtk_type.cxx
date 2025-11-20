#include <cstring>

#include <mdcraft/io/vtk/vtk_type.h>

namespace mdcraft::io {

VtkType::VtkType(const char* value) {
    if (!std::strcmp(value, "Undefined")) {
        m_value = Undefined;
    } else if (!std::strcmp(value, "Int8")) {
        m_value = Int8;
    } else if (!std::strcmp(value, "Int16")) {
        m_value = Int16;
    } else if (!std::strcmp(value, "Int32")) {
        m_value = Int32;
    } else if (!std::strcmp(value, "Int64")) {
        m_value = Int64;
    } else if (!std::strcmp(value, "UInt8")) {
        m_value = UInt8;
    } else if (!std::strcmp(value, "UInt16")) {
        m_value = UInt16;
    } else if (!std::strcmp(value, "UInt32")) {
        m_value = UInt32;
    } else if (!std::strcmp(value, "UInt64")) {
        m_value = UInt64;
    } else if (!std::strcmp(value, "Float32")) {
        m_value = Float32;
    } else if (!std::strcmp(value, "Float64")) {
        m_value = Float64;
    } else {
        m_value = Undefined;
    }
}

VtkType::VtkType(const std::string& value)
    : VtkType(value.c_str()) { }

const char* VtkType::to_str() const {
    switch (m_value) {
        case Int8:    return "Int8";
        case Int16:   return "Int16";
        case Int32:   return "Int32";
        case Int64:   return "Int64";

        case UInt8:   return "UInt8";
        case UInt16:  return "UInt16";
        case UInt32:  return "UInt32";
        case UInt64:  return "UInt64";

        case Float32: return "Float32";
        case Float64: return "Float64";

        default: return "Undefined";
    }
}

size_t VtkType::size() const {
    switch (m_value) {
        case VtkType::Int64:
        case VtkType::UInt64:
        case VtkType::Float64:
            return 8;
        case VtkType::Int32:
        case VtkType::UInt32:
        case VtkType::Float32:
            return 4;
        case VtkType::Int16:
        case VtkType::UInt16:
            return 2;
        case VtkType::Int8:
        case VtkType::UInt8:
            return 1;
        default:
            return 0;
    }
}

} // namespace mdcraft::io