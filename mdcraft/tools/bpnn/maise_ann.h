#pragma once

#include <mdcraft/tools/bpnn/base.h>
#include <mdcraft/tools/parser/maise_model.h>

namespace mdcraft {
namespace bpnn {

using ::mdcraft::parser::MaiseModelParser;

class MaiseBPNN : public IBPNN {
public:
	MaiseBPNN(std::string_view filename);

	const std::vector<std::vector<double>>&         vRmin() const;
	const std::vector<std::vector<double>>&         vRmax() const;
	const std::vector<std::vector<double>>&         range() const;

	~MaiseBPNN() {}

private:
	std::unique_ptr<MaiseModelParser> 				parser;
};

} // bpnn
} // mdcraft
