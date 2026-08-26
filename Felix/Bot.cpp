#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <dpp/dpp.h>
#include "Felix.h"
#include "bmp/stb_image_write.h"
#include "bmp/gif.h"

struct CommandReplyVisitor
{
	const dpp::message_create_t& event;

	void operator()(std::monostate)
	{
		dpp::message m;
		m.content = "----> done";
		event.reply(m.set_flags(dpp::m_ephemeral));
	}

	void operator()(const Felix::ErrorMessage& message)
	{
		dpp::message m;
		m.content = "--X-> failed: " + message.string;
		event.reply(m.set_flags(dpp::m_ephemeral));
	}

	void operator()(const std::string& res)
	{
		dpp::embed embed;
		embed.add_field("", res);
		dpp::message m{ event.msg.channel_id, embed };
		event.reply(m.set_flags(dpp::m_ephemeral));
	}

	void operator()(const Felix::ImageList& images)
	{
		if (images.empty()) {
			dpp::message m;
			m.content = "--X-> failed: no images";
			event.reply(m.set_flags(dpp::m_ephemeral));
			return;
		}

		// png
		if (images.size() == 1) {
			auto& image = images[0];

			std::vector<uint8_t> rgb;
			rgb.resize(image.resolution * image.resolution * 3);
			for (size_t y = 0; y != image.resolution; ++y) {
				for (size_t x = 0; x != image.resolution; ++x) {
					size_t ry = image.resolution - y - 1;
					size_t offset = (y * image.resolution + x);
					auto& pixel = image[offset];
					if (pixel.size() != 3) {
						continue;
					}
					size_t roffset = (ry * image.resolution + x);
					size_t rgb_offset = roffset * 3;
					rgb[rgb_offset + 0] = pixel[0] * 255;
					rgb[rgb_offset + 1] = pixel[1] * 255;
					rgb[rgb_offset + 2] = pixel[2] * 255;
				}
			}

			int size = 0;
			auto png = stbi_write_png_to_mem(
				rgb.data(),
				image.resolution * 3,
				image.resolution,
				image.resolution,
				3,
				&size);
			std::string_view content{ (char*)png, (size_t)size };

			dpp::embed embed;
			embed.set_image("attachment://image.png");
			dpp::message m{ event.msg.channel_id, embed };
			m.add_file("image.png", content);
			free(png);

			event.reply(m);
			return;
		}

		// gif
		{
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

			dpp::embed embed;
			embed.set_image("attachment://graph.gif");
			dpp::message m{ event.msg.channel_id, embed };
			std::string_view content{ (const char*)data.data(), data.size() };
			m.add_file("graph.gif", content);
			event.reply(m);
			return;
		}
		
	};

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
	std::string BOT_TOKEN;
	{
		const char* token_file_path = getenv("FELIXBOT_TOKEN_FILE");
		if (!token_file_path) token_file_path = "UNSET";
		std::ifstream token_file{ token_file_path };
		token_file >> BOT_TOKEN;
	}

	auto intents = dpp::i_default_intents | dpp::i_message_content;
	dpp::cluster bot(BOT_TOKEN, intents );

	bot.on_log(dpp::utility::cout_logger());

	bot.on_slashcommand([](const dpp::slashcommand_t& event) {
		if (event.command.get_command_name() == "ping") {
			event.reply("Pong!");
		}
	});

	Felix f;
	bot.on_message_create([&](const dpp::message_create_t& event) {
		if (event.msg.author.is_bot()) return;
		if (event.msg.content.starts_with(">>")) {
			std::string expr = event.msg.content.substr(2);
#if 0
			if (expr.starts_with("help")) {
				expr = expr.substr(5);
				dpp::embed embed;
				/**/ if (expr.starts_with("functions")) {
					embed.set_title("Math functions");
					std::string
						power_functions_str,
						trig_functions_str,
						for_functions_str,
						value_functions_str,
						logic_functions_str,
						other_functions_str,
						experimental_functions_str;

					power_functions_str += "```\n";
					power_functions_str += "     exp(x) -> exponent\n";
					power_functions_str += "      ln(x) -> natural logarithm\n";
					power_functions_str += "    sqrt(x) -> square root\n";
					power_functions_str += "inv_sqrt(x) -> inverse square root\n";
					power_functions_str += " root[y](x) -> y-th root of x\n";
					power_functions_str += "  log[y](x) -> logarithm of x to base y\n";
					power_functions_str += "```";
					embed.add_field("Power functions", power_functions_str, false);

					trig_functions_str += "```\n";
					trig_functions_str += "    cos(x)       cot(x)\n";
					trig_functions_str += "   cosh(x)      coth(x)\n";
					trig_functions_str += " arccos(x)    arccot(x)\n";
					trig_functions_str += "arccosh(x)   arccoth(x)\n";
					trig_functions_str += "    sin(x)       tan(x)\n";
					trig_functions_str += "   sinh(x)      tanh(x)\n";
					trig_functions_str += " arcsin(x)    arctan(x)\n";
					trig_functions_str += "arcsinh(x)   arctanh(x)\n";
					trig_functions_str += "   sin1(x)\n";
					trig_functions_str += "   cos1(x)\n";
					trig_functions_str += "```";
					embed.add_field("Trigonometric functions", trig_functions_str, false);

					for_functions_str += "```\n";
					for_functions_str += "S{t;begin;end}[f(t)] -> Sum\n";
					for_functions_str += "       FD{t;x}[f(t)] -> Forward Difference\n";
					for_functions_str += "     FD{t;x;n}[f(t)] -> n-th Forward Difference\n";
					for_functions_str += "       BD{t;x}[f(t)] -> Backward Difference\n";
					for_functions_str += "     BD{t;x;n}[f(t)] -> n-th Backward Difference\n";
					for_functions_str += "P{t;begin;end}[f(t)] -> Product\n";
					for_functions_str += "      R{t;x;n}[f(t)] -> Return\n";
					for_functions_str += "      I{t;a;b}[f(t)] -> Integral\n";
					for_functions_str += "        D{t;x}[f(t)] -> Derivative\n";
					for_functions_str += "      D{t;x;n}[f(t)] -> n-th Derivative\n";
					for_functions_str += "   Iexp{t;a;b}[f(t)] -> Integral along exp\n";
					for_functions_str += "```";
					embed.add_field("FOR-like functions", for_functions_str, false);

					value_functions_str += "```\n";
					value_functions_str += "    abs(x) -> absolute value\n";
					value_functions_str += "inv_abs(x) -> inverse absolute value\n";
					value_functions_str += "    arg(x) -> argument of x\n";
					value_functions_str += "   sign(x) -> normalized number\n";
					value_functions_str += "     Re(x) -> real part\n";
					value_functions_str += "     Im(x) -> imaginary part\n";
					value_functions_str += "  floor(x)\n";
					value_functions_str += "   ceil(x)\n";
					value_functions_str += "  round(x)\n";
					value_functions_str += "   norm(x) | normalize(x)\n";
					value_functions_str += "```";
					embed.add_field("Value functions", value_functions_str, false);

					logic_functions_str += "```\n";
					logic_functions_str += "  exist(x) -> if number exist returns 0 else 1\n";
					logic_functions_str += "   grid(x) -> if number on the grid returns 0 else (0 < ...)\n";
					logic_functions_str += "if{x}[a;b] -> if x is 0 then returns a else returns b\n";
					logic_functions_str += "```";
					embed.add_field("Logic functions", logic_functions_str, false);

					other_functions_str += "```\n";
					other_functions_str += "  gamma(x) -> gamma function\n";
					other_functions_str += " fctI(x,n) -> n-th integral of factorial\n";
					other_functions_str += "rand | rand(...)\n";
					other_functions_str += "```";
					embed.add_field("Other functions", other_functions_str, false);

					experimental_functions_str += "```\n";
					experimental_functions_str += "Poly{t;f(t);x}[a;b;...] -> Polynomial\n";
					experimental_functions_str += "                 zbf(x) -> zeta(x)/x!\n";
					experimental_functions_str += "             USumN(x,n) -> Undefined sum of x^n\n";
					experimental_functions_str += "```";
					embed.add_field("Experimental functions", experimental_functions_str, false);
				}
				else if (expr.starts_with("commands")) {
					embed.set_title("Felix commands");
				}
				else if (expr.starts_with("operators")) {
					embed.set_title("Math operators");
				}
				else {
					embed.set_title("Help menu");
				}
				dpp::message m{ event.msg.channel_id, embed };
				event.reply(m);
				return;
			}
#endif
			CommandReplyVisitor crv(event);
			std::visit(crv, f.run_command(expr));
			return;
		}
		if (event.msg.content.starts_with(">")) {
			std::string expr = event.msg.content.substr(1);
			dpp::embed embed;
			embed.add_field(expr, f.evaluate(expr));
			dpp::message m{ event.msg.channel_id, embed };
			event.reply(m);
			return;
		}
	});

	bot.on_ready([&bot](const dpp::ready_t& event) {
		if (dpp::run_once<struct register_bot_commands>()) {
			auto delete_cb = bot.co_global_bulk_command_delete();
			delete_cb.sync_wait();

			bot.global_command_create(dpp::slashcommand("ping", "Ping pong!", bot.me.id));
		}
	});

	bot.start(dpp::st_wait);
}
