//	v1.1.0

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"

bool run_command(std::string c) {
	return false;
}

int main() {
	std::string expression;
	for (;;) {
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			if (!run_command(expression.substr(1)))
				std::cout << "unknown command";
		}
		else {
			Function f("f", {}, parse_expr(expression));
			math::complex res = f.return_value();
			std::cout << " -> " << res.R << "  " << res.i << "i";
		}
		std::cout << "\n\n";
	}
}