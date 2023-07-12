//	v1.2.3

#include <iostream>
#include "include.h"
#include "object/object.h"
#include "bmp/BMPWriter.h"

bool run_command(std::string c) {
	std::vector<std::string> args;
	for (;;) {
		int index = c.find(' ');
		if (-1 < index) {
			args.push_back(c.substr(0, index));
			c = c.substr(index + 1);
		}
		else {
			args.push_back(c);
			break;
		}
	}

	if (args[0] == "add") {
		Function f;

		int bracket_index = args[1].find('(');
		if (-1 < bracket_index) {
			f.name = args[1].substr(0, bracket_index);
			std::string arguments = args[1].substr(bracket_index + 1);
			while (arguments.length() > 0) {
				int index1 = arguments.find(','),
					index2 = arguments.find(')');
				/**/ if (-1 < index1) f.args.push_back(Variable(arguments.substr(0, index1), 0));
				else if (-1 < index2) f.args.push_back(Variable(arguments.substr(0, index2), 0));
				else /*------------*/ f.args.push_back(Variable(arguments, 0));
				arguments = arguments.substr(f.args[f.args.size() - 1].name.length() + 1);
			}
		}
		else f.name = args[1];

		std::string expr = "";
		for (int i = 2; i < args.size(); i++) expr += args[i];
		f.objects = parse_expr(expr);

		functions_list.push_back(f);
		return true;
	}
	if (args[0] == "remove") {
		for (int i = 0; i < functions_list.size(); i++)
			if (functions_list[i].name == args[1])
				functions_list.erase(functions_list.begin() + i);
		return true;
	}
	if (args[0] == "list") {
		for (auto f : functions_list) {
			std::cout << " -> " << f.name;
			if (f.args.size() > 0) {
				std::cout << "(";
				for (auto a : f.args) std::cout << a.name << ",";
				std::cout << "\b)";
			}
			std::cout << "\n";
		}
		return true;
	}

	return false;
}

int main() {
	std::string expression;
	for (;;) {
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			if (run_command(expression.substr(1)))	std::cout << " ----> done";
			else									std::cout << " --x-> unknown command";
		}
		else {
			Function f("f", {}, parse_expr(expression));
			math::complex res = f.return_value();
			std::cout << " -> " << math::toString(res);
		}
		std::cout << "\n\n";
	}
}