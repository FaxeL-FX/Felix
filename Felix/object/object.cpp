#include "object.h"
#include <chrono>

std::vector<Function> functions_list;

mobj::mathObjs nameToType(std::string token) {
	if (token == "_Const") return mobj::mathObjs::_Const;

	if (token == "+") return mobj::mathObjs::_add;
	if (token == "*") return mobj::mathObjs::_mul;
	if (token == "/") return mobj::mathObjs::_div;
	if (token == "%") return mobj::mathObjs::_mod;
	if (token == "^") return mobj::mathObjs::_pow;
	
	if (token == "abs") return mobj::mathObjs::abs;
	if (token == "inv_abs") return mobj::mathObjs::inv_abs;
	if (token == "arg") return mobj::mathObjs::arg;
	if (token == "sign") return mobj::mathObjs::sign;
	if (token == "Re") return mobj::mathObjs::Re;
	if (token == "Im") return mobj::mathObjs::Im;
	if (token == "exist") return mobj::mathObjs::exist;
	if (token == "floor") return mobj::mathObjs::floor;
	if (token == "round") return mobj::mathObjs::round;
	if (token == "grid") return mobj::mathObjs::grid;
	
	if (token == "exp") return mobj::mathObjs::exp;
	if (token == "ln") return mobj::mathObjs::ln;
	if (token == "sqrt") return mobj::mathObjs::sqrt;
	if (token == "inv_sqrt") return mobj::mathObjs::inv_sqrt;
	if (token == "root") return mobj::mathObjs::root;
	if (token == "log") return mobj::mathObjs::log;
	
	if (token == "cos") return mobj::mathObjs::cos;
	if (token == "cosh") return mobj::mathObjs::cosh;
	if (token == "arccos") return mobj::mathObjs::arccos;
	if (token == "arccosh") return mobj::mathObjs::arccosh;
	if (token == "sin") return mobj::mathObjs::sin;
	if (token == "sinh") return mobj::mathObjs::sinh;
	if (token == "arcsin") return mobj::mathObjs::arcsin;
	if (token == "arcsinh") return mobj::mathObjs::arcsinh;
	if (token == "cot" || token == "ctg") return mobj::mathObjs::cot;
	if (token == "coth" || token == "ctgh") return mobj::mathObjs::coth;
	if (token == "arccot" || token == "arcctg") return mobj::mathObjs::arccot;
	if (token == "arccoth" || token == "arcctgh") return mobj::mathObjs::arccoth;
	if (token == "tan" || token == "tg") return mobj::mathObjs::tan;
	if (token == "tanh" || token == "tgh") return mobj::mathObjs::tanh;
	if (token == "arctan" || token == "arctg") return mobj::mathObjs::arctan;
	if (token == "arctanh" || token == "arctgh") return mobj::mathObjs::arctanh;

	if (token == "gamma") return mobj::mathObjs::gamma;
	if (token == "fctIntegral" || token == "fctI") return mobj::mathObjs::fctIntegral;

	if (token == "S" || token == "Sum") return mobj::mathObjs::Sum;
	if (token == "P" || token == "Product") return mobj::mathObjs::Product;
	if (token == "R" || token == "Return") return mobj::mathObjs::Return;
	if (token == "I" || token == "Integral") return mobj::mathObjs::Integral;
	if (token == "D" || token == "Derivative") return mobj::mathObjs::Derivative;
	if (token == "Iexp" || token == "IntegralAlongExp") return mobj::mathObjs::IntegralAlongExp;

	return mobj::mathObjs::_Default;
}

std::vector<Object> parse_expr(std::string expr) {
	if (expr.length() == 0) return { Object(-1, -1, "_Const", {}, 0) };
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
	parse_obj(&objects, tokens, 0, tokens.size() - 1);
	for (int i = 0; i < objects.size(); i++) {
		for (auto fnc : functions_list)
			if (objects[i].name == fnc.name) {
				objects[i].fncID = fnc.id;
				break;
			}
	}

	return objects;
}
std::string parse_token(std::string expr) {
	if (expr.length() == 0) return " ";
	int type = 0;	//	types: 0->any , 1->number , 2->word
	char c = std::tolower(expr[0]);
	if ('0' <= expr[0] && expr[0] <= '9' || expr[0] == '.') type = 1;
	else if ('a' <= c && c <= 'z') type = 2;

	if (type == 0) switch (expr[0]) {
	case('+', '-', '*', '/', '%', '^', '!', '(', ')', '[', ']', '{', '}', '=', ',', ';'):
		return expr.substr(0, 1);
	case('<', '>'):
		if ((1 < expr.length()) && expr[1] == '=')	return expr.substr(0, 2);
		else										return expr.substr(0, 1);
	}
	int i = 1;
	for (i; i < expr.length(); i++) {
		c = std::tolower(expr[i]);
		if (('0' <= c && c <= '9' || c == '.') && type == 1 ||
			('a' <= c && c <= 'z' || c == '_') && type == 2) continue;
		break;
	}
	return expr.substr(0, i);
}
void parse_obj(std::vector<Object> *objects, std::vector<std::string> tokens, int begin, int end) {
	objects->push_back(Object());
	std::string token = tokens[end];
	int index = end, objIndex = objects->size() - 1;
	if (begin < 0 || tokens.size() <= end) {
		(*objects)[objIndex].type = mobj::mathObjs::_Error;
		return;
	}
	(*objects)[objIndex].type = mobj::mathObjs::_Default;

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
		(*objects)[objIndex].type = nameToType(token);
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, index + 1, end);
		return;
	}
	if (token == "-") {
		if (index == begin) {
			(*objects)[objIndex].type = mobj::mathObjs::_uMinus;
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		else {
			(*objects)[objIndex].type = mobj::mathObjs::_dif;
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, begin, index - 1);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		return;
	}
	if (token == "!") {
		(*objects)[objIndex].type = mobj::mathObjs::_fct;
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		return;
	}
	if ('0' <= token[0] && token[0] <= '9' || token[0] == '.') {
		(*objects)[objIndex].type = mobj::mathObjs::_Const;
		(*objects)[objIndex].value = std::strtold(token.c_str(), 0);
		return;
	}

	/**/ if (token == "pi") {
		(*objects)[objIndex].type = mobj::mathObjs::_Const;
		(*objects)[objIndex].value = math::pi;
		return;
	}
	else if (token == "e") {
		(*objects)[objIndex].type = mobj::mathObjs::_Const;
		(*objects)[objIndex].value = math::exp(1);
		return;
	}
	else if (token == "i") {
		(*objects)[objIndex].type = mobj::mathObjs::_Const;
		(*objects)[objIndex].value = math::i;
		return;
	}
	if (index + 1 < tokens.size()) {
		if (tokens[index + 1] == "(") {
			(*objects)[objIndex].name = token;
			(*objects)[objIndex].type = nameToType(token);

			int start_i = index + 2;
			index = brackets(tokens, index + 1, true);
			for (int i = start_i; i <= index; i++) {
				if (tokens[i] == ")" || tokens[i] == ",") {
					(*objects)[objIndex].arg_indexes.push_back(objects->size());
					parse_obj(objects, tokens, start_i, i - 1);
					start_i = i + 1;
					continue;
				}
				if (tokens[i] == "(" || tokens[i] == "[" || tokens[i] == "{")
					i = brackets(tokens, i, true);
			}
			return;
		}
		if (tokens[index + 1] == "[" && tokens[brackets(tokens, index + 1, true) + 1] == "(") {
			(*objects)[objIndex].name = token;
			(*objects)[objIndex].type = nameToType(token);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));

			index = brackets(tokens, index + 1, true);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
			return;
		}
		if (tokens[index + 1] == "{" && tokens[brackets(tokens, index + 1, true) + 1] == "[") {
			(*objects)[objIndex].name = token;
			(*objects)[objIndex].type = nameToType(token);

			int start_i = index + 2;
			index = brackets(tokens, index + 1, true);
			for (int i = start_i; i <= index; i++) {
				if (tokens[i] == "}" || tokens[i] == ";") {
					(*objects)[objIndex].arg_indexes.push_back(objects->size());
					parse_obj(objects, tokens, start_i, i - 1);
					start_i = i + 1;
					continue;
				}
				if (tokens[i] == "(" || tokens[i] == "[" || tokens[i] == "{")
					i = brackets(tokens, i, true);
			}

			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
			return;
		}
	}
	(*objects)[objIndex].name = token;

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
	//	(...) | ...! | ...;... | ...,...
	if (token1 == "(" || token2 == ")" || token2 == "!" || token1 == ";" || token2 == ";" || token1 == "," || token2 == ",") return false;
	//	fnc()
	if (('a' <= token1[0] && token1[0] <= 'z' || 'A' <= token1[0] && token1[0] <= 'Z' || token1[0] == '_') && token2 == "(") return false;

	return true;
}

math::complex Object::return_value(std::vector<Object>* objects, std::vector<Variable>* args) {
	std::vector<math::complex> args_results;
	for (auto i : this->arg_indexes) args_results.push_back((*objects)[i].return_value(objects, args));

	switch (arg_indexes.size()) {
	case(0): switch (this->type) {
		case(mobj::mathObjs::_Error): return math::exp(3000);
		case(mobj::mathObjs::_Const): return this->value;
	}
	case(1): switch (this->type) {
		case(mobj::mathObjs::_uMinus): return -args_results[0];
		case(mobj::mathObjs::_fct): return math::fct(args_results[0]);

		case(mobj::mathObjs::abs): return math::abs(args_results[0]);
		case(mobj::mathObjs::inv_abs): return math::inv_abs(args_results[0]);
		case(mobj::mathObjs::arg): return math::arg(args_results[0]);
		case(mobj::mathObjs::sign): return math::normalize(args_results[0]);
		case(mobj::mathObjs::Re): return ((math::complex_linear)args_results[0]).R;
		case(mobj::mathObjs::Im): return ((math::complex_linear)args_results[0]).i;
		case(mobj::mathObjs::exist): return !math::exist(args_results[0]);
		case(mobj::mathObjs::floor): return math::floor(args_results[0]);
		case(mobj::mathObjs::round): return math::floor(args_results[0] + math::complex(0.5, 0.5));
		case(mobj::mathObjs::grid): {
			math::complex_linear res = args_results[0] - math::floor(args_results[0] + math::complex(0.5, 0.5));
			if (res.R < 0) res.R = -res.R;
			if (res.i < 0) res.i = -res.i;
			if (res.R > res.i) return res.i;
			return res.R;
		}

		case(mobj::mathObjs::exp): return math::exp(args_results[0]);
		case(mobj::mathObjs::ln): return math::ln(args_results[0]);
		case(mobj::mathObjs::sqrt): return math::sqrt(args_results[0]);
		case(mobj::mathObjs::inv_sqrt): return math::inv_sqrt(args_results[0]);

		case(mobj::mathObjs::cos): return math::cos(args_results[0]);
		case(mobj::mathObjs::cosh): return math::cosh(args_results[0]);
		case(mobj::mathObjs::arccos): return math::arccos(args_results[0]);
		case(mobj::mathObjs::arccosh): return math::arccosh(args_results[0]);
		case(mobj::mathObjs::sin): return math::sin(args_results[0]);
		case(mobj::mathObjs::sinh): return math::sinh(args_results[0]);
		case(mobj::mathObjs::arcsin): return math::arcsin(args_results[0]);
		case(mobj::mathObjs::arcsinh): return math::arcsinh(args_results[0]);
		case(mobj::mathObjs::cot): return math::cot(args_results[0]);
		case(mobj::mathObjs::coth): return math::coth(args_results[0]);
		case(mobj::mathObjs::arccot): return math::arccot(args_results[0]);
		case(mobj::mathObjs::arccoth): return math::arccoth(args_results[0]);
		case(mobj::mathObjs::tan): return math::tan(args_results[0]);
		case(mobj::mathObjs::tanh): return math::tanh(args_results[0]);
		case(mobj::mathObjs::arctan): return math::arctan(args_results[0]);
		case(mobj::mathObjs::arctanh): return math::arctanh(args_results[0]);

		case(mobj::mathObjs::gamma): return math::fct(args_results[0] - 1);
		case(mobj::mathObjs::fctIntegral): return math::fctIntegral(args_results[0], 0);
	}
	case(2): switch (this->type) {
		case(mobj::mathObjs::_add): return args_results[0] + args_results[1];
		case(mobj::mathObjs::_dif): return args_results[0] - args_results[1];
		case(mobj::mathObjs::_mul): return args_results[0] * args_results[1];
		case(mobj::mathObjs::_div): return args_results[0] / args_results[1];
		case(mobj::mathObjs::_pow): return math::pow(args_results[0], args_results[1]);
		case(mobj::mathObjs::_mod): return args_results[0] % args_results[1];

		case(mobj::mathObjs::root): return math::pow(args_results[1], 1 / args_results[0]);
		case(mobj::mathObjs::log): return math::ln(args_results[1]) / math::ln(args_results[0]);
	}
	case(3): switch (this->type) {
		case(mobj::mathObjs::Derivative): {
			const int n = 256;
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1] + 0.5 / n));
			math::complex res = (*objects)[this->arg_indexes[2]].return_value(objects, args);
			(*args)[var_index].value = args_results[1] - 0.5 / n;
			res = res - (*objects)[this->arg_indexes[2]].return_value(objects, args);
			args->erase(args->begin() + var_index);
			return res * n;
		}
	}
	case(4): switch (this->type) {
		case(mobj::mathObjs::Sum): {
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::complex
				difference = args_results[2] - args_results[1],
				step = math::normalize(difference),
				res = 0;
			int repeats = math::abs(difference) + 1;
			for (int i = 0; i < repeats; i++) {
				res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args);
				(*args)[var_index].value = (*args)[var_index].value + step;
			}
			args->erase(args->begin() + var_index);
			return res;
		}
		case(mobj::mathObjs::Product): {
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::complex
				difference = args_results[2] - args_results[1],
				step = math::normalize(difference),
				res = 1;
			int repeats = math::abs(difference) + 1;
			for (int i = 0; i < repeats; i++) {
				res = res * (*objects)[this->arg_indexes[3]].return_value(objects, args);
				(*args)[var_index].value = (*args)[var_index].value + step;
			}
			args->erase(args->begin() + var_index);
			return res;
		}
		case(mobj::mathObjs::Return): {
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1]));
			int repeats = ((math::complex_linear)args_results[2]).R;
			for (int i = 0; i < repeats; i++) {
				(*args)[var_index].value = (*objects)[this->arg_indexes[3]].return_value(objects, args);
			}
			math::complex res = (*args)[var_index].value;
			args->erase(args->begin() + var_index);
			return res;
		}

		case(mobj::mathObjs::Integral): {
			const int n = 256;
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::complex res = 0, dx = (args_results[2] - args_results[1]) / n;
			for (int k = 0; k < n; k++) {
				res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args);
				(*args)[var_index].value = (*args)[var_index].value + dx;
			}
			args->erase(args->begin() + var_index);
			return res * dx;
		}
		case(mobj::mathObjs::IntegralAlongExp): {
			const int n = 256;
			int var_index = args->size();
			args->push_back(Variable(arg_indexes.size(), (*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::complex res = 0, dx = (math::ln(args_results[2] / args_results[1])) / n, x = math::ln(args_results[1]), dt;
			for (int k = 0; k < n; k++) {
				dt = math::exp(x + dx) - math::exp(x);
				res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * dt;
				x = x + dx;
				(*args)[var_index].value = math::exp(x);
			}
			args->erase(args->begin() + var_index);
			return res;
		}
	}
	}

	if (this->type == mobj::mathObjs::rand) {
		if (args_results.size() == 0) return math::rand(std::chrono::steady_clock::now().time_since_epoch().count() % 4096);
		return math::rand(std::chrono::steady_clock::now().time_since_epoch().count() % 4096, args_results);
	}
	for (auto a : (*args))
		if (a.id == this->argID) return a.value;
	for (auto f : functions_list)
		if (f.id == this->fncID) {
			for (int i = 0; i < f.args.size() && i < args_results.size(); i++)
				f.args[i].value = args_results[i];
			return f.return_value();
		}
	return 0;
}

math::complex value(std::vector<Object> objects, std::vector<Variable> args) {
	return objects[0].return_value(&objects, &args);
}