#pragma once
#include "include.h"
#include "object/object.h"

struct Felix {
	math::complex plot_center = 0;
	math::real plot_radius = 8;

	std::vector<Function> functions_list;


	bool run_command(std::string command);
};