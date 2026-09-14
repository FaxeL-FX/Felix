#pragma once
#include <memory>

#include "include.h"
#include "object/object.h"
#include "text/text_image.h"

struct Felix {
	struct ErrorMessage {
		std::string string;
	};

	struct Image : std::vector<std::array<float, 3>> {
		Image(const std::vector<std::array<float, 3>> &data, int resolution)
			: std::vector<std::array<float, 3>>{data},
			  resolution{resolution}
		{
		}
		int resolution;
	};

	struct ImageList : std::vector<Image> {
		uint32_t gif_delay;
	};

	struct PrintTask
	{
		bool do_work();
		bool done();
		ImageList images();

		~PrintTask();
		struct Impl;
		std::unique_ptr<Impl> I;
	};

	using CommandResult = std::variant<
		std::monostate,
		ErrorMessage,
		std::string,
		ImageList
	>;

	math::complex plot_center;
	long double plot_radius = 8.0;

	std::vector<Function> functions_list;


	CommandResult run_command(std::string command);
	std::string evaluate(std::string expression);
};
