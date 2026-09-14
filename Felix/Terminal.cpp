#include "Felix.h"
#include <Windows.h>
#include <fstream>
#include <thread>
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

	void operator()(Felix::PrintTask &print_task)
	{
		auto worker = [&print_task](){ while (!print_task.done()) print_task.do_work();};
		std::vector<std::thread> thread_pool;
		size_t prc = std::thread::hardware_concurrency();
		if (prc == 0) prc++;
		for (size_t i = 0; i != prc; ++i) thread_pool.emplace_back(worker);

		std::cout << " render [" << std::to_string(100.0).substr(1) << "%]\b\b";
		while (!print_task.done()) {
			float percent = 100 + 100 * print_task.progress();
			std::cout << "\b\b\b\b\b\b\b\b\b" << std::to_string(percent).substr(1);
			std::this_thread::sleep_for(std::chrono::milliseconds{100});
		}
		for (auto &T : thread_pool) T.join();
		std::cout << "\b\b\b\b\b\b\b\b\b\bcomplete     \n";

		auto images = print_task.images();
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