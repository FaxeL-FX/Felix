#include "fixed_point.h"
#include <print>


int main() {
	math::fixed_point<32> 
		a{63.0},
		b{63};

	std::println("a = {}", a.toString(20));
	std::println("b = {}", b.toString(20));

	std::println("a+b = {}", (a + b).toString(20));
	std::println("a-b = {}", (a - b).toString(20));
	std::println("a*b = {}", (a * b).toString(20));
	std::println("a/b = {}", (a / b).toString(20));

	std::println("exp(a) = {}", math::exp(a).toString(20));
	std::println("exp(b) = {}", math::exp(b).toString(20));
}