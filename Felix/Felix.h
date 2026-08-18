#pragma once
#include "include.h"
#include "object/object.h"

struct Felix {
	struct ErrorMessage {
		std::string string;
	};

	struct Image : std::vector<std::vector<float>> {
		Image(const std::vector<std::vector<float>> &data, int resolution)
			: std::vector<std::vector<float>>{data},
			  resolution{resolution}
		{
		}
		int resolution;
	};

	struct ImageList : std::vector<Image> {
		uint32_t gif_delay;
	};

	using CommandResult = std::variant<
		std::monostate,
		ErrorMessage,
		std::string,
		ImageList
	>;

	math::complex plot_center = 0;
	math::real plot_radius = 8;

	std::vector<Function> functions_list;


	CommandResult run_command(std::string command);
	std::string evaluate(std::string expression);
};
