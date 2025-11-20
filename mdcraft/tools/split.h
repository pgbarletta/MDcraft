#pragma once

#include <vector>
#include <string>
#include <regex>

namespace mdcraft {
namespace tools {

inline std::vector<std::string>
custom_split_regex(const std::string& str, const std::string& delim)
{
	std::vector<std::string> out;
	std::regex rgx(delim);
	std::sregex_token_iterator it(str.begin(), str.end(), rgx, -1), end;
	if (*it == "") ++it;
	for (; it != end; ++it) out.push_back(*it);
	return out;
}

inline std::vector<std::string>
custom_split(const std::string& str, const std::string& delim)
{
	if (delim.empty())
		throw std::runtime_error("A null-delimiter is geven\n");
	if (delim.find(" ") != std::string::npos)
		throw std::runtime_error("A delimiter must not contain white spaces\n");
	return custom_split_regex(str, delim);
}

inline std::vector<std::string>
custom_split(const std::string& str, const char delim=' ')
{
	if (delim == ' ')
		return custom_split_regex(str, "\\s+");
	return custom_split_regex(str, std::string(1, delim));
}

} // tools
} // mdcraft