
#include "Felix.h"
#include <Windows.h>

int main() {
	Felix felix;

	SetConsoleOutputCP(CP_UTF8);
	std::string expression;
	for (;;) {
		std::cout << " ";
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			if (felix.run_command(expression.substr(1)))	std::cout << " ----> done";
			else											std::cout << " --X-> failed";
		}
		else {
			auto objs = parse_expr(expression);
			std::vector<Variable> vars = {};
			math::number answer = value(objs, vars, felix.functions_list);
			std::cout << " -> " << answer.toString();
		}
		std::cout << "\n\n";
	}
}