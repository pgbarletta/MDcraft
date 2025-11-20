#pragma once
#ifdef DEBUG

#include <vector>
#include <iostream>

namespace mdcraft {
namespace tools {

template <typename T>
void print(T v) {
	std::cout << v << std::endl;
}

template <typename T>
void print(const std::vector<T>& v) {
	for (auto&& e : v) std::cout << e << " ";
	std::cout << std::endl;
}

} // tools
} // mdcraft

#endif