
#include "Felix.h"
#include <Windows.h>
#include <fstream>
#include <iterator>
#include "bmp/BMPWriter.h"
#include "bmp/gif.h"

#include "include.h"

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

		if (images.size() == 1) {
			std::ofstream image_file{"graph.bmp"};
			auto &first_image = images[0];
			std::string bmp_image =
				BMPWriter::make_bmp_image(first_image, first_image.resolution);
			image_file << bmp_image;
			std::cout << " ----> done (to graph.bmp)";
			return;
		}

		auto resolution = images[0].resolution;

		std::vector<uint8_t> data;
		GifWriterCallback cb;
		cb.udata = &data;
		cb.write = [](void *udata, void *data, size_t size) {
			std::vector<uint8_t> &vec = *static_cast<std::vector<uint8_t>*>(udata);
			uint8_t *start = (uint8_t *)data;
			uint8_t *end = start + size;
			std::copy(start, end, std::back_inserter(vec));
		};

		GifWriter writer = {};
		GifBegin(&writer, cb, resolution, resolution, images.gif_delay);
		for (auto &image : images) {
			write_gif_frame(&writer, image, resolution, images.gif_delay);
		}
		GifEnd(&writer);
		{
		    std::ofstream of{"graph.gif", std::ios::binary};
		    of.write((char*)data.data(), data.size());
		}
		std::cout << " ----> done (to graph.gif)";
	}

	void write_gif_frame(
		GifWriter *writer,
		const std::vector<std::vector<float>> &img,
		int resolution,
		uint32_t delay)
	{
		std::vector<uint8_t> image;
		image.resize(resolution * resolution * 4);
		for (int yy = 0; yy < resolution; ++yy) {
			for (int xx = 0; xx < resolution; ++xx) {
				size_t offset_image = ((resolution - yy - 1) * resolution + xx);
				size_t offset = (yy * resolution + xx);
				auto &r = image[offset_image * 4 + 0];
				auto &g = image[offset_image * 4 + 1];
				auto &b = image[offset_image * 4 + 2];
				auto &a = image[offset_image * 4 + 3];
				r = img[offset][0] * 255;
				g = img[offset][1] * 255;
				b = img[offset][2] * 255;
				a = 255;
			}
		}
		GifWriteFrame(writer, image.data(), resolution, resolution, delay);
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