#include <dpp/dpp.h>
#include <print>
#include "Felix.h"
#include "ImageListBinary.hh"

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

		ImageListBinary bin{images};
		std::string attach_name = std::format("attachment://{}", bin.file_name);
		dpp::embed embed;
		embed.set_image(attach_name);
		dpp::message m{event.msg.channel_id, embed};
		m.add_file(bin.file_name, std::string_view{(char*)bin.data(),bin.size()});
		event.reply(m);
	};
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
#if 1
			if (expr.starts_with("help")) {
				expr = expr.substr(4);
				dpp::embed embed;
				/**/ if (expr.starts_with(" functions")) {
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
				else if (expr.starts_with(" commands")) {
					embed.set_title("Felix commands");
				}
				else if (expr.starts_with(" operators")) {
					embed.set_title("Math operators");
				}
				else {
					embed.set_title("Help menu");
					embed.add_field("", ">>help functions\n>>help commands\n>>help operators");
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
