#pragma once

#include <fstream>
#include <string_view>

namespace mdcraft {
namespace parser {

class IFileContent;

class IParser {
public:
	IParser(std::string_view filename);

	virtual ~IParser();

	virtual void readHead();
    virtual void readBody();

    virtual bool headRead() const;
    virtual bool bodyRead() const;

protected:
	std::ifstream inputFile_{""};

	bool is_head_read_{ false };
	bool is_body_read_{ false };
};

class FileOpenException : public std::runtime_error {
public:
    FileOpenException(std::string_view filename);

    std::string_view getFilename() const;

    const char* what() const noexcept override;
    // TODO: add method what() to fix segfault

private:
    std::string_view filename_;
};

class IFileContent {

};

} // parser
} // mdcraft
