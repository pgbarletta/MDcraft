#pragma once

#include <vector>
#include <memory>
#include <string_view>

#include <mdcraft/tools/bpnn/base.h>
#include <mdcraft/tools/parser/maise_model.h>

namespace mdcraft {
namespace bpnn {

using ::mdcraft::parser::MaiseModelParser;

class MaiseBPDescriptor : public IBPDescriptor {
public:
	MaiseBPDescriptor(std::string_view filename);

	~MaiseBPDescriptor();

private:
	std::unique_ptr<MaiseModelParser> 				parser;

	void fillShiftTable();
};

} // bpnn
} // mdcraft
