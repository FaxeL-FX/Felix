#include "Felix.h"
#include <Windows.h>
#include <fstream>
#include "ImageListBinary.hh"

struct CommandResultVisitor
{
	void operator()(std::monostate)
	{
		std::cout << " ----> done";
	}

	void operator()(const Felix::ErrorMessage &message)
	{
		std::cout << " --X-> failed: " << message.string;
	}

	void operator()(const std::string &res)
	{
		std::cout << res;
	}

	void operator()(const Felix::ImageList &images)
	{
		if (images.empty()) {
			std::cout << " ----> Empty image list";
			return;
		}

		ImageListBinary bin{images};
		std::ofstream image_file{bin.file_name, std::ios::binary};
		image_file.write((char*)bin.data(), bin.size());
		std::cout << std::format(" ----> done (to {})", bin.file_name);
	}
};

int main() {
	Felix felix;

	SetConsoleOutputCP(CP_UTF8);
	std::string expression;
	for (;;) {
		std::cout << " ";
		std::getline(std::cin, expression);
		if (expression[0] == '>') {
			Felix::CommandResult res = felix.run_command(expression.substr(1));
			std::visit(CommandResultVisitor{}, res);
		}
		else {
			std::cout << " -> " << felix.evaluate(expression);
		}
		std::cout << "\n\n";
	}
}