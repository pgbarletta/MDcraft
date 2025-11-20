#include <mdcraft/tools/parser/base.h>

namespace mdcraft {
namespace parser {

FileOpenException::FileOpenException(std::string_view filename)
  : filename_(filename)
  , std::runtime_error(std::string("Failed to open file: ") + filename_.data()) {}

std::string_view FileOpenException::getFilename() const { return filename_; }

const char* FileOpenException::what() const noexcept {
  return (std::string("Failed to open file: ") + filename_.data()).data();
}

IParser::IParser(std::string_view filename)
 : inputFile_(filename.data(), std::ios_base::in)
{
	if (!inputFile_.is_open()) {
		inputFile_.close();
		throw FileOpenException(filename);
	}
}

void IParser::readHead() { return; }
void IParser::readBody() { return; }

bool IParser::headRead() const { return is_head_read_; }
bool IParser::bodyRead() const { return is_body_read_; }

IParser::~IParser() {
	inputFile_.close();
}

} // parser
} // mdcraft