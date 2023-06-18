#include "object.h"
#include <chrono>

std::vector<Function> functions_list;

std::vector<Object> parse_expr(std::string expr) {
	if (expr.length() == 0) return { Object("_Const", {}, 0) };
	std::vector<std::string> tokens;
	while (expr.length() > 0) {
		std::string token = parse_token(expr);
		if (token != " ") tokens.push_back(token);
		expr = expr.substr(token.length());
	}
	for (int i = 1; i < tokens.size(); i++)
		if (need_mul(tokens[i - 1], tokens[i]))
			tokens.insert(tokens.begin() + i, "*");

	std::vector<Object> objects;
	std::vector<Object>* objs = &objects;
	parse_obj(objs, tokens, 0, tokens.size() - 1);
	return objects;
}
std::string parse_token(std::string expr) {
	if (expr.length() == 0) return " ";
	int type = 0;	//	types: 0->any , 1->number , 2->word
	char c = std::tolower(expr[0]);
	if ('0' <= expr[0] && expr[0] <= '9' || expr[0] == '.' || expr[0] == ',') type = 1;
	else if ('a' <= c && c <= 'z') type = 2;

	if (type == 0) switch (expr[0]) {
	case('+', '-', '*', '/', '%', '^', '!', '(', ')', '[', ']', '{', '}', '='):
		return expr.substr(0, 1);
	case('<', '>'):
		if ((1 < expr.length()) && expr[1] == '=')	return expr.substr(0, 2);
		else										return expr.substr(0, 1);
	}
	int i = 1;
	for (i; i < expr.length(); i++) {
		c = std::tolower(expr[i]);
		if (('0' <= c && c <= '9' || c == '.' || c == ',') && type == 1 ||
			('a' <= c && c <= 'z' || c == '_') && type == 2) continue;
		break;
	}
	return expr.substr(0, i);
}
void parse_obj(std::vector<Object> *objects, std::vector<std::string> tokens, int begin, int end) {
	if (begin < 0) begin = 0;
	if (tokens.size() <= end) end = tokens.size() - 1;
	objects->push_back(Object());
	std::string token = tokens[end];
	int index = end, objIndex = objects->size() - 1;
	(*objects)[objIndex].name = "_Undefined";

	for (int i = end; begin <= i && i <= end; i--) {
		if (tokens[index] == "(" || tokens[index] == "[" || tokens[index] == "{" || tokens[i] == ")" || tokens[i] == "]" || tokens[i] == "}") {
			if (i == end && brackets(tokens, i, false) == begin) {
				objects->pop_back();
				return parse_obj(objects, tokens, begin + 1, end - 1);
			}
			else i = brackets(tokens, i, false);
			continue;
		}
		if (priority(tokens[i]) <= priority(token)) continue;
		token = tokens[i];
		index = i;
	}
	if (token == "+" || token == "*" || token == "/" || token == "%" || token == "^") {
		(*objects)[objIndex].name = token;
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, index + 1, end);
		return;
	}
	if (token == "-") {
		(*objects)[objIndex].name = token;
		if (index == begin) {
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		else {
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, begin, index - 1);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		return;
	}
	if (token == "!") {
		(*objects)[objIndex].name = token;
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		return;
	}
	if ('0' <= token[0] && token[0] <= '9' || token[0] == '.' || token[0] == ',') {
		(*objects)[objIndex].name = "_Const";
		(*objects)[objIndex].value = std::strtold(token.c_str(), 0);	//	need to change for "flex_float" instead "long double"
		return;
	}

	/**/ if (token == "pi") {
		(*objects)[objIndex].name = "_Const";
		(*objects)[objIndex].value = math::pi;
	}
	else if (token == "e") {
		(*objects)[objIndex].name = "_Const";
		(*objects)[objIndex].value = math::exp(1);
	}
	else if (token == "i") {
		(*objects)[objIndex].name = "_Const";
		(*objects)[objIndex].value = math::i;
	}
	else if (token == "rand") {
		(*objects)[objIndex].name = "_Rand";
	}
	else {
		std::vector<std::string>
			fnc_1 = {
				"floor",
				"exp",
				"ln",
				"sqrt",
				"inv_sqrt",
				"cosh",
				"sinh",
				"coth", "ctgh",
				"tanh", "tgh",
				"arccosh",
				"arcsinh",
				"arccoth", "arcctgh",
				"arctanh", "arctgh",
				"cos",
				"sin",
				"cot", "ctg",
				"tan", "tg",
				"arccos",
				"arcsin",
				"arccot", "arcctg",
				"arctan", "arctg",
				"rand",
			}, 
			fnc_2 = {
				"root",
				"log",
			},
			fnc_3 = {
				"S", "Sum",
				"P", "Product",
				"R", "Return",
			};
		int fnc = -1;
		for (auto s : fnc_1) if (token == s) fnc = 1;
		for (auto s : fnc_2) if (token == s) fnc = 2;
		for (auto s : fnc_3) if (token == s) fnc = 3;
		switch (fnc) {
		case(1):
			if (tokens[index + 1] == "(") {
				(*objects)[objIndex].name = token;
				(*objects)[objIndex].arg_indexes.push_back(objects->size());
				parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
				return;
			}
			return;
		case(2):
			if (tokens[index + 1] == "[" && tokens[brackets(tokens, index + 1, true) + 1] == "(") {
				(*objects)[objIndex].name = token;
				(*objects)[objIndex].arg_indexes.push_back(objects->size());
				parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));

				index = brackets(tokens, index + 1, true);
				(*objects)[objIndex].arg_indexes.push_back(objects->size());
				parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
				return;
			}
			return;
		case(3):
			if (tokens[index + 1] == "{" && tokens[brackets(tokens, index + 1, true) + 1] == "(") {
				return;

				(*objects)[objIndex].name = token;
				//(*objects)[objIndex].arg_indexes.push_back(objects->size());
				//parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));

				index = brackets(tokens, index + 1, true);
				(*objects)[objIndex].arg_indexes.push_back(objects->size());
				parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
				return;
			}
			return;
		default:
			(*objects)[objIndex].name = token;
		}
	}
}

int brackets(std::vector<std::string> tokens, int bracket_index, bool forward) {
	char antibracket = antibr(tokens[bracket_index][0]);
	if (forward) {
		for (int index = bracket_index + 1; 0 <= index && index < tokens.size(); index++) {
			if (tokens[index][0] == antibracket) return index;
			if (tokens[index] == "(" || tokens[index] == "[" || tokens[index] == "{")
				index = brackets(tokens, index, forward);
		}
		return tokens.size() - 1;
	}
	for (int index = bracket_index - 1; 0 <= index && index < tokens.size(); index--) {
		if (tokens[index][0] == antibracket) return index;
		if (tokens[index] == ")" || tokens[index] == "]" || tokens[index] == "}")
			index = brackets(tokens, index, forward);
	}
	return 0;
}
char antibr(char b) {
	switch (b) {
	case('('): return ')';
	case('['): return ']';
	case('{'): return '}';
	case(')'): return '(';
	case(']'): return '[';
	case('}'): return '{';
	default: return b;
	}
}
int priority(std::string token) {
	if (token == "+" || token == "-") return 5;
	else if (token == "*" || token == "/" || token == "%") return 4;
	else if (token == "^") return 3;
	else if (token == "!") return 2;
	else if (token == "(" || token == "[" || token == "{" || token == ")" || token == "]" || token == "}") return 0;
	return 1;
}
int find_token(std::vector<std::string> tokens, int token_index, bool forward) {
	if (token_index < 1 || tokens.size() <= token_index) return 0;
	if (forward) {
		for (int index = token_index + 1; index < tokens.size(); index++) {
			if (tokens[index] == tokens[token_index]) return index;
			switch (tokens[index][0]) {
			case('(', '[', '{'): index = brackets(tokens, index, forward);
			}
		}
		return tokens.size() - 1;
	}
	for (int index = token_index - 1; 0 <= index; index--) {
		if (tokens[index] == tokens[token_index]) return index;
		switch (tokens[index][0]) {
		case(')', ']', '}'): index = brackets(tokens, index, forward);
		}
	}
	return 0;
}
bool need_mul(std::string token1, std::string token2) {
	//	+ - * / ^
	if (priority(token1) > 2 || priority(token2) > 2) return false;
	//	[...] | {...}
	if (token1 == "[" || token1 == "]" || token1 == "{" || token1 == "}" || token2 == "[" || token2 == "]" || token2 == "{" || token2 == "}") return false;
	//	(...) | ...!
	if (token1 == "(" || token2 == ")" || token2 == "!") return false;
	//	fnc()
	if (('a' <= token1[0] && token1[0] <= 'z' || token1[0] == '_') && token2 == "(") return false;

	return true;
}

math::complex Object::return_value(std::vector<Object>* objects, std::vector<Variable>* args) {
	std::vector<math::complex> args_results;
	for (auto i : this->arg_indexes) args_results.push_back((*objects)[i].return_value(objects, args));

	switch (args_results.size()) {
	case(0): {
		if (this->name == "_Const") return this->value;
		if (this->name == "_Rand") return math::rand(std::chrono::steady_clock::now().time_since_epoch().count() % 4096);
	}
	case(1): {
		if (this->name == "-") return -args_results[0];
		if (this->name == "!") return math::fct(args_results[0]);

		if (this->name == "floor") return math::floor(args_results[0]);

		if (this->name == "exp") return math::exp(args_results[0]);
		if (this->name == "ln") return math::ln(args_results[0]);
		if (this->name == "sqrt") return math::sqrt(args_results[0]);
		if (this->name == "inv_sqrt") return math::inv_sqrt(args_results[0]);

		if (this->name == "cos") return math::cos(args_results[0]);
		if (this->name == "cosh") return math::cosh(args_results[0]);
		if (this->name == "arccos") return math::arccos(args_results[0]);
		if (this->name == "arccosh") return math::arccosh(args_results[0]);
		if (this->name == "sin") return math::sin(args_results[0]);
		if (this->name == "sinh") return math::sinh(args_results[0]);
		if (this->name == "arcsin") return math::arcsin(args_results[0]);
		if (this->name == "arcsinh") return math::arcsinh(args_results[0]);
		if (this->name == "cot" || this->name == "ctg") return math::cot(args_results[0]);
		if (this->name == "coth" || this->name == "ctgh") return math::coth(args_results[0]);
		if (this->name == "arccot" || this->name == "arcctg") return math::arccot(args_results[0]);
		if (this->name == "arccoth" || this->name == "arcctgh") return math::arccoth(args_results[0]);
		if (this->name == "tan" || this->name == "tg") return math::tan(args_results[0]);
		if (this->name == "tanh" || this->name == "tgh") return math::tanh(args_results[0]);
		if (this->name == "arctan" || this->name == "arctg") return math::arctan(args_results[0]);
		if (this->name == "arctanh" || this->name == "arctgh") return math::arctanh(args_results[0]);
	}
	case(2): {
		if (this->name == "+") return args_results[0] + args_results[1];
		if (this->name == "-") return args_results[0] - args_results[1];
		if (this->name == "*") return args_results[0] * args_results[1];
		if (this->name == "/") return args_results[0] / args_results[1];
		if (this->name == "%") return args_results[0] % args_results[1];
		if (this->name == "^") return math::pow(args_results[0], args_results[1]);

		if (this->name == "root") return math::pow(args_results[1], 1 / args_results[0]);
		if (this->name == "log") return math::ln(args_results[1]) / math::ln(args_results[0]);
	}
	default: return 0;
	}
}
