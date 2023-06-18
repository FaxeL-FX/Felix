//	v1.0.2

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"

int main() {
	std::string expression;
	for (;;) {
		std::getline(std::cin, expression);
		Function f("f", {}, parse_expr(expression));
		math::complex res = f.return_value();
		std::cout << " -> " << res.R << "  " << res.i << "i\n\n";
	}
}