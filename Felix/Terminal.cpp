#include "Felix.h"
#include <Windows.h>
#include <algorithm>
#include <fstream>
#include <thread>
#include "ImageListBinary.hh"

struct PrintTaskInProgress : std::optional<Felix::PrintTask>
{
	void set(Felix::PrintTask &&new_task)
	{
		if (has_value()) abort_task();
		emplace(std::move(new_task));
		aborted = false;
		manager_thread = std::jthread(&PrintTaskInProgress::manager_routine, this);
	}

	void abort_task()
	{
		aborted = true;
		if (manager_thread.joinable()) manager_thread.join();
	}

	void manager_routine()
	{
		Felix::PrintTask& task = value();
		aborted = false;

		std::vector<std::jthread> workers;
		size_t n_workers = std::thread::hardware_concurrency();
		n_workers = std::max<size_t>(1, n_workers);
		for (size_t i = 0; i != n_workers; ++i) workers.emplace_back(&PrintTaskInProgress::worker_routine, this);

		std::string report_buf;
		while (not aborted and not task.done()) {
			auto report = std::format(" render [{:.6f}%]", task.progress() * 100);
			if (report.size() > report_buf.size()) report_buf.resize(report.size());
			std::ranges::fill(report_buf, ' ');
			std::copy(report.begin(), report.end(), report_buf.begin());
			std::cout << '\r' << report_buf << '\r';
			std::cout.flush();
			std::this_thread::sleep_for(std::chrono::milliseconds{ 10 });
		}

		workers.clear(); // jthread::join

		std::ranges::fill(report_buf, ' ');
		std::cout << '\r' << report_buf << '\r';

		if (aborted) {
			std::cout << " render aborted\n ";
			return;
		}

		auto images = task.images();
		if (images.empty()) {
			std::cout << " ----> Empty image list\n ";
			return;
		}

		ImageListBinary bin{ images };
		std::ofstream image_file{ bin.file_name, std::ios::binary };
		image_file.write((char*)bin.data(), bin.size());
		std::cout << std::format(" ----> done (to {})\n ", bin.file_name);
	}

	void worker_routine()
	{
		Felix::PrintTask &task = value();
		while (not aborted and not task.done()) task.do_work();
	}

	std::atomic<bool> aborted = false;
	std::jthread manager_thread;
};

struct CommandResultVisitor
{
	PrintTaskInProgress &active_print_task;

	void operator()(std::monostate)
	{
		std::cout << " ----> done\n";
	}

	void operator()(const Felix::ErrorMessage &message)
	{
		std::cout << " --X-> failed: " << message.string << "\n";
	}

	void operator()(const std::string &res)
	{
		std::cout << res << "\n";
	}

	void operator()(Felix::PrintTask &print_task)
	{
		active_print_task.set(std::move(print_task));
	}
};

int main() {
	Felix felix;
	PrintTaskInProgress active_print_task;

	SetConsoleOutputCP(CP_UTF8);
	std::string expression;
	for (;;) {
		std::cout << " ";
		std::getline(std::cin, expression);
		active_print_task.abort_task();
		if (expression[0] == '>') {
			Felix::CommandResult res = felix.run_command(expression.substr(1));
			std::visit(CommandResultVisitor{active_print_task}, res);
		}
		else {
			std::cout << " -> " << felix.evaluate(expression) << "\n";
		}
	}
}