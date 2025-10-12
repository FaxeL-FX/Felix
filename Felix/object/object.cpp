#include "object.h"
#include <chrono>

std::vector<Function> functions_list;

ObjType nameToType(std::string token) {
	if ('0' <= token[0] && token[0] <= '9' || token[0] == '.') return ObjType::_Const;

	if (token == "+") return ObjType::_add;
	if (token == "*") return ObjType::_mul;
	if (token == "/") return ObjType::_div;
	if (token == "%") return ObjType::_mod;
	if (token == "^") return ObjType::_pow;

	if (token == "!" || token == "fct") return ObjType::_fct;
	
	if (token == "exp")			return ObjType::_exp;
	if (token == "ln")			return ObjType::_ln;
	if (token == "sqrt")		return ObjType::_sqrt;
	if (token == "inv_sqrt")	return ObjType::_inv_sqrt;
	if (token == "root")		return ObjType::_root;
	if (token == "log")			return ObjType::_log;
	
	if (token == "cos")		return ObjType::_cos;
	if (token == "cosh")	return ObjType::_cosh;
	if (token == "arccos")	return ObjType::_arccos;
	if (token == "arccosh") return ObjType::_arccosh;
	if (token == "sin")		return ObjType::_sin;
	if (token == "sinh")	return ObjType::_sinh;
	if (token == "arcsin")	return ObjType::_arcsin;
	if (token == "arcsinh") return ObjType::_arcsinh;
	if (token == "cot" || token == "ctg")			return ObjType::_cot;
	if (token == "coth" || token == "ctgh")			return ObjType::_coth;
	if (token == "arccot" || token == "arcctg")		return ObjType::_arccot;
	if (token == "arccoth" || token == "arcctgh")	return ObjType::_arccoth;
	if (token == "tan" || token == "tg")			return ObjType::_tan;
	if (token == "tanh" || token == "tgh")			return ObjType::_tanh;
	if (token == "arctan" || token == "arctg")		return ObjType::_arctan;
	if (token == "arctanh" || token == "arctgh")	return ObjType::_arctanh;

	if (token == "sin1")	return ObjType::_sin1;
	if (token == "cos1")	return ObjType::_cos1;
	if (token == "Binom")	return ObjType::_Binom;

	if (token == "S" || token == "Sum")						return ObjType::_Sum;
	if (token == "FDiff" || token == "ForwardDifference")	return ObjType::_ForwardDifference;
	if (token == "P" || token == "Product")					return ObjType::_Product;
	if (token == "R" || token == "Return")					return ObjType::_Return;
	if (token == "I" || token == "Integral")				return ObjType::_Integral;
	if (token == "D" || token == "Derivative")				return ObjType::_Derivative;
	if (token == "Iexp" || token == "IntegralAlongExp")		return ObjType::_IntegralAlongExp;
	if (token == "Poly" || token == "Polynomial")			return ObjType::_Polynomial;

	if (token == "abs")		return ObjType::_abs;
	if (token == "inv_abs") return ObjType::_inv_abs;
	if (token == "arg")		return ObjType::_arg;
	if (token == "sign")	return ObjType::_sign;
	if (token == "Re")		return ObjType::_Re;
	if (token == "Im")		return ObjType::_Im;
	if (token == "floor")	return ObjType::_floor;
	if (token == "ceil")	return ObjType::_ceil;
	if (token == "round")	return ObjType::_round;

	if (token == "exist")	return ObjType::_exist;
	if (token == "grid")	return ObjType::_grid;

	if (token == "gamma")							return ObjType::_gamma;
	if (token == "fctIntegral" || token == "fctI")	return ObjType::_fctIntegral;
	if (token == "H")								return ObjType::_Harmonic;
	if (token == "zeta")							return ObjType::_zeta;


	if (token == "rand")	return ObjType::_rand;

	return ObjType::_Default;
}

std::vector<Object> parse_expr(std::string expr) {
	if (expr.length() == 0) return { Object("_Const", {}, 0) };
	std::vector<std::string> tokens;
	while (expr.length() > 0) {
		std::string token = parse_token(expr);
		if (token[0] != ' ') tokens.push_back(token);
		expr = expr.substr(token.length());
	}
	//for (auto t : tokens) std::cout << t << "\n";

	for (int i = 1; i < tokens.size(); i++)
		if (need_mul(tokens[i - 1], tokens[i]))
			tokens.insert(tokens.begin() + i, "*");

	std::vector<Object> objects;
	parse_obj(&objects, tokens, 0, tokens.size() - 1);
	for (int i = 0; i < objects.size(); i++) {
		objects[i].id = 0;
		for (char c : objects[i].name)
			objects[i].id = (objects[i].id << 1) + (int)c;
	}

	return objects;
}
std::string parse_token(std::string expr) {
	if (expr.length() == 0) return " ";
	int type = 0, i = 1;	//	types: 0->any , 1->number , 2->word
	char c = std::tolower(expr[0]);
	if ('0' <= expr[0] && expr[0] <= '9' || expr[0] == '.') type = 1;
	else if ('a' <= c && c <= 'z' || c == '_') type = 2;

	// number
	if (type == 1) {
		for (i; i < expr.length(); i++) {
			if (('0' <= expr[i] && expr[i] <= '9' || expr[i] == '.')) continue;
			break;
		}
	}

	// word
	if (type == 2) {
		for (i; i < expr.length(); i++) {
			c = std::tolower(expr[i]);
			if ('a' <= c && c <= 'z' || c == '_' || '0' <= c && c <= '9') continue;
			break;
		}
	}
	return expr.substr(0, i);
}
void parse_obj(std::vector<Object> *objects, std::vector<std::string> tokens, int begin, int end) {
	if (begin < 0 || tokens.size() <= end) return;

	int index = end, objIndex = objects->size();
	std::string token = tokens[end];

	if (tokens[end] == ")" || tokens[end] == "]" || tokens[end] == "}")
		if (brackets(tokens, end, false) == begin)
			return parse_obj(objects, tokens, begin + 1, end - 1);

	for (int i = end; begin <= i && i <= end; i--) {
		if (tokens[i] == ")" || tokens[i] == "]" || tokens[i] == "}") {
			i = brackets(tokens, i, false);
			continue;
		}
		if (priority(tokens[i]) <= priority(token)) continue;
		index = i;
		token = tokens[index];
	}

	objects->push_back(Object());
	(*objects)[objIndex].type = ObjType::_Default;

	if (token == "+" || token == "*" || token == "/" || token == "%" || token == "^") {
		if (index <= begin || end <= index) {
			(*objects)[objIndex].type = ObjType::_Error;
			return;
		}
		(*objects)[objIndex].type = nameToType(token);
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, index + 1, end);
		return;
	}
	if (token == "-") {
		if (index == begin) {
			if (end <= index) {
				(*objects)[objIndex].type = ObjType::_Error;
				return;
			}
			(*objects)[objIndex].type = ObjType::_uMinus;
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		else {
			if (index < begin || end <= index) {
				(*objects)[objIndex].type = ObjType::_Error;
				return;
			}
			(*objects)[objIndex].type = ObjType::_dif;
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, begin, index - 1);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, end);
		}
		return;
	}
	if (token == "!") {
		if (index <= begin) {
			(*objects)[objIndex].type = ObjType::_Error;
			return;
		}
		(*objects)[objIndex].type = ObjType::_fct;
		(*objects)[objIndex].arg_indexes.push_back(objects->size());
		parse_obj(objects, tokens, begin, index - 1);
		return;
	}
	if ('0' <= token[0] && token[0] <= '9' || token[0] == '.') {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = std::strtold(token.c_str(), 0);
		return;
	}

	/**/ if (token == "pi") {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = math::pi;
		return;
	}
	else if (token == "e") {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = math::exp(1);
		return;
	}
	else if (token == "i") {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = math::i;
		return;
	}
#if infsimIsHere
	else if (token == "inf") {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = math::infinity;
		return;
	}
	else if (token == "o") {
		(*objects)[objIndex].type = ObjType::_Const;
		(*objects)[objIndex].value = math::zero;
		return;
	}
#endif

	(*objects)[objIndex].name = token;
	(*objects)[objIndex].type = nameToType(token);
	if (index + 1 < tokens.size()) {
		if (tokens[index + 1] == "(") {
			int start_i = index + 2;
			index = brackets(tokens, index + 1, true);
			for (int i = start_i; i <= index; i++) {
				if (tokens[i] == ")" || tokens[i] == ";" || tokens[i] == ",") {
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
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));

			index = brackets(tokens, index + 1, true);
			(*objects)[objIndex].arg_indexes.push_back(objects->size());
			parse_obj(objects, tokens, index + 1, brackets(tokens, index + 1, true));
			return;
		}
		if (tokens[index + 1] == "{" && tokens[brackets(tokens, index + 1, true) + 1] == "[") {
			int start_i = index + 2;
			index = brackets(tokens, index + 1, true);
			for (int i = start_i; i <= index; i++) {
				if (tokens[i] == "}" || tokens[i] == ";" || tokens[i] == ",") {
					(*objects)[objIndex].arg_indexes.push_back(objects->size());
					parse_obj(objects, tokens, start_i, i - 1);
					start_i = i + 1;
					continue;
				}
				if (tokens[i] == "(" || tokens[i] == "[" || tokens[i] == "{")
					i = brackets(tokens, i, true);
			}
			start_i = index + 2;
			index = brackets(tokens, index + 1, true);
			for (int i = start_i; i <= index; i++) {
				if (tokens[i] == "]" || tokens[i] == ";" || tokens[i] == ",") {
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
	}
}

int brackets(std::vector<std::string> tokens, int bracket_index, bool forward) {
	char antibracket = antibr(tokens[bracket_index][0]);
	if (forward)
		for (int index = bracket_index + 1; 0 <= index && index < tokens.size(); index++) {
			if (tokens[index][0] == antibracket) return index;
			if (tokens[index] == "(" || tokens[index] == "[" || tokens[index] == "{")
				index = brackets(tokens, index, forward);
		}
	else
		for (int index = bracket_index - 1; 0 <= index && index < tokens.size(); index--) {
			if (tokens[index][0] == antibracket) return index;
			if (tokens[index] == ")" || tokens[index] == "]" || tokens[index] == "}")
				index = brackets(tokens, index, forward);
		}
	return bracket_index;
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
bool need_mul(std::string token1, std::string token2) {
	//	+ - * / ^
	if (priority(token1) > 2 || priority(token2) > 2) return false;

	if (token1 == "(" || token1 == "[" || token1 == "{" || token1 == ";" || token1 == ",") return false;
	if (token2 == ")" || token2 == "]" || token2 == "}" || token2 == ";" || token2 == "," || token2 == "!") return false;
	// fnc[]() | fnc{}[]
	if (token1 == "]" && token2 == "(" || token1 == "}" && token2 == "[") return false;
	//	fnc()
	token1[0] = std::tolower(token1[0]);
	if (('a' <= token1[0] && token1[0] <= 'z' || token1[0] == '_') && (token2 == "(" || token2 == "[" || token2 == "{")) return false;

	return true;
}

math::number Object::return_value(std::vector<Object>* objects, std::vector<Variable>* args) {
	std::vector<math::number> args_results;
	for (auto i : this->arg_indexes) args_results.push_back((*objects)[i].return_value(objects, args));

	switch (arg_indexes.size()) {
	case(0): switch (this->type) {
		case(ObjType::_Error): return math::exp(3000);
		case(ObjType::_Const): return this->value;
	}
	case(1): switch (this->type) {
		case(ObjType::_uMinus): return -args_results[0];
		case(ObjType::_fct): return math::fct(args_results[0]);

		case(ObjType::_exp): return math::exp(args_results[0]);
		case(ObjType::_ln): return math::ln(args_results[0]);
		case(ObjType::_sqrt): return math::sqrt(args_results[0]);
		case(ObjType::_inv_sqrt): return math::inv_sqrt(args_results[0]);

		case(ObjType::_cos): return math::cos(args_results[0]);
		case(ObjType::_cosh): return math::cosh(args_results[0]);
		case(ObjType::_arccos): return math::arccos(args_results[0]);
		case(ObjType::_arccosh): return math::arccosh(args_results[0]);
		case(ObjType::_sin): return math::sin(args_results[0]);
		case(ObjType::_sinh): return math::sinh(args_results[0]);
		case(ObjType::_arcsin): return math::arcsin(args_results[0]);
		case(ObjType::_arcsinh): return math::arcsinh(args_results[0]);
		case(ObjType::_cot): return math::cot(args_results[0]);
		case(ObjType::_coth): return math::coth(args_results[0]);
		case(ObjType::_arccot): return math::arccot(args_results[0]);
		case(ObjType::_arccoth): return math::arccoth(args_results[0]);
		case(ObjType::_tan): return math::tan(args_results[0]);
		case(ObjType::_tanh): return math::tanh(args_results[0]);
		case(ObjType::_arctan): return math::arctan(args_results[0]);
		case(ObjType::_arctanh): return math::arctanh(args_results[0]);
		case(ObjType::_sin1): return math::sin1(args_results[0]);
		case(ObjType::_cos1): return math::cos1(args_results[0]);

		case(ObjType::_abs): return math::abs(args_results[0]);
		case(ObjType::_inv_abs): return math::inv_abs(args_results[0]);
		case(ObjType::_arg): return math::arg(args_results[0]);
		case(ObjType::_sign): return math::normalize(args_results[0]);
		case(ObjType::_Re): return Re(args_results[0]);
		case(ObjType::_Im): return Im(args_results[0]);

		case(ObjType::_floor): return math::floor(args_results[0]);
		case(ObjType::_ceil): return math::ceil(args_results[0]);
		case(ObjType::_round): return math::floor(args_results[0] + (math::number)(math::complex(0.5, 0.5)));

		case(ObjType::_exist): return !math::exist(args_results[0]);
		case(ObjType::_grid): return math::grid(args_results[0]);

		case(ObjType::_gamma): return math::fct(args_results[0] - 1);
		case(ObjType::_Harmonic): return math::Harmonic(args_results[0]);
		case(ObjType::_zeta): return math::zeta(args_results[0]);
	}
	case(2): switch (this->type) {
		case(ObjType::_add): return args_results[0] + args_results[1];
		case(ObjType::_dif): return args_results[0] - args_results[1];
		case(ObjType::_mul): return args_results[0] * args_results[1];
		case(ObjType::_div): return args_results[0] / args_results[1];
		case(ObjType::_pow): return math::pow(args_results[0], args_results[1]);
		case(ObjType::_mod): return args_results[0] % args_results[1];

		case(ObjType::_root): return math::pow(args_results[1], 1 / args_results[0]);
		case(ObjType::_log): return math::ln(args_results[1]) / math::ln(args_results[0]);

		case(ObjType::_fctIntegral): return math::fctIntegral(args_results[0], args_results[1]);

		case(ObjType::_Binom): return math::Binom(args_results[0], args_results[1]);
	}
	case(3): switch (this->type) {
		case(ObjType::_Derivative): {
			const int n = 256;
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1] + 0.5 / n));
			math::number res = (*objects)[this->arg_indexes[2]].return_value(objects, args);
			(*args)[var_index].value = args_results[1] - 0.5 / n;
			res = res - (*objects)[this->arg_indexes[2]].return_value(objects, args);
			args->erase(args->begin() + var_index);
			return res * n;
		}
		case(ObjType::_ForwardDifference): {
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::number res = (*objects)[this->arg_indexes[2]].return_value(objects, args);
			(*args)[var_index].value = args_results[1] - 1;
			res = res - (*objects)[this->arg_indexes[2]].return_value(objects, args);
			args->erase(args->begin() + var_index);
			return res;
		}
	}
	case(4): switch (this->type) {
		case(ObjType::_Sum):
		case(ObjType::_Product): {
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::number
				res = (*objects)[this->arg_indexes[3]].return_value(objects, args),
				difference = args_results[2] - args_results[1];
#if infsimIsHere
			long long
				repeats = math::floor(difference.getNum(math::acch).R);
			math::number
				diffPart = math::number(difference.getNum(math::acch).R - repeats, difference.getNum(math::acch).i),
				step = math::normalize(diffPart.getNum(math::acch));
#else
			long long
				repeats = math::floor(difference.R);
			math::number
				diffPart = math::number(difference.R - repeats, difference.i),
				step = math::normalize(diffPart);
#endif
			long double
				radius = math::abs(diffPart),
				radFract = radius - math::floor(radius);
			math::number
				angle1 = math::mul_i((math::number)math::arccos(0.5 * radFract)),
				angle2 = math::mul_i((math::number)math::arccos(0.5 * (1 - radFract))),
				step1 =  step * math::exp( angle1),
				step2 =  step * math::exp(-angle1),
				step3 = -step * math::exp( angle2),
				step4 = -step * math::exp(-angle2);

			(*args)[var_index].value = (*args)[var_index].value + 0.5 * (1 - step);
			/**/ if (this->type == ObjType::_Sum) {
				for (int i = 0; i < (long long)radius; i++) {
					(*args)[var_index].value = (*args)[var_index].value + step;
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step;
				}
				if (radFract != 0) {
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step + step1);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step1 * 0.5 * (1 - radFract);
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step1 + step2);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step2 * 0.5 * (1 - radFract);
					(*args)[var_index].value = (*args)[var_index].value - step1;
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step2 * 0.5 * (1 - radFract);
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step2 + step1);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step1 * 0.5 * (1 - radFract);

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step - step1) - step2;
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step * radFract;

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step + step3);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step3 * 0.5 * radFract;
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step3 + step4);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step4 * 0.5 * radFract;
					(*args)[var_index].value = (*args)[var_index].value - step3;
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step4 * 0.5 * radFract;
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step4 + step3);
					res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * step3 * 0.5 * radFract;

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step3 - 1);
				}
				if (repeats < 0) {
					for (int i = 0; repeats < i; i--) {
						res = res - (*objects)[this->arg_indexes[3]].return_value(objects, args);
						(*args)[var_index].value = (*args)[var_index].value - 1;
					}
				}
				else {
					for (int i = 0; i < repeats; i++) {
						(*args)[var_index].value = (*args)[var_index].value + 1;
						res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args);
					}
				}
			}
			else if (this->type == ObjType::_Product) {
				for (int i = 0; i < (int)radius; i++) {
					(*args)[var_index].value = (*args)[var_index].value + step;
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step);
				}
				if (radFract != 0) {
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step + step1);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step1 * 0.5 * (1 - radFract));
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step1 + step2);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step2 * 0.5 * (1 - radFract));
					(*args)[var_index].value = (*args)[var_index].value - step1;
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step2 * 0.5 * (1 - radFract));
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step2 + step1);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step1 * 0.5 * (1 - radFract));

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step - step1) - step2;
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step * radFract);

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step + step3);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step3 * 0.5 * radFract);
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step3 + step4);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step4 * 0.5 * radFract);
					(*args)[var_index].value = (*args)[var_index].value - step3;
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step4 * 0.5 * radFract);
					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step4 + step3);
					res = res * math::pow((*objects)[this->arg_indexes[3]].return_value(objects, args), step3 * 0.5 * radFract);

					(*args)[var_index].value = (*args)[var_index].value + 0.5 * (step3 - 1);
				}
				if (repeats < 0) {
					for (int i = 0; repeats < i; i--) {
						res = res / (*objects)[this->arg_indexes[3]].return_value(objects, args);
						(*args)[var_index].value = (*args)[var_index].value - 1;
					}
				}
				else {
					for (int i = 0; i < repeats; i++) {
						(*args)[var_index].value = (*args)[var_index].value + 1;
						res = res * (*objects)[this->arg_indexes[3]].return_value(objects, args);
					}
				}
			}

			args->erase(args->begin() + var_index);
			return res;
		}
		case(ObjType::_Return): {
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			int repeats = math::round(((math::complex)args_results[2]).R);
			for (int i = 0; i < repeats; i++) {
				(*args)[var_index].value = (*objects)[this->arg_indexes[3]].return_value(objects, args);
			}
			math::number res = (*args)[var_index].value;
			args->erase(args->begin() + var_index);
			return res;
		}

		case(ObjType::_Integral): {
			const int n = 1024;
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::number res = 0, dx = (args_results[2] - args_results[1]) / n;
			for (int k = 0; k < n; k++) {
				res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args);
				(*args)[var_index].value = (*args)[var_index].value + dx;
			}
			args->erase(args->begin() + var_index);
			return res * dx;
		}
		case(ObjType::_IntegralAlongExp): {
			const int n = 1024;
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::number res = 0, dx = (math::ln(args_results[2] / args_results[1])) / n, x = math::ln(args_results[1]), dt;
			for (int k = 0; k < n; k++) {
				dt = math::exp(x + dx) - math::exp(x);
				res = res + (*objects)[this->arg_indexes[3]].return_value(objects, args) * dt;
				x = x + dx;
				(*args)[var_index].value = math::exp(x);
			}
			args->erase(args->begin() + var_index);
			return res;
		}

		case(ObjType::_ForwardDifference): {
			int var_index = args->size();
			args->push_back(Variable((*objects)[this->arg_indexes[0]].name, args_results[1]));
			math::number res = (*objects)[this->arg_indexes[3]].return_value(objects, args);

			int n = ((math::complex)(*objects)[this->arg_indexes[2]].return_value(objects, args)).R;
			for (int k = 1; k <= n; k++) {
				(*args)[var_index].value = (*args)[var_index].value - 1;
				res = res + (1 - 2 * (k % 2)) * 
					math::factorial(n) / (math::factorial(k) * math::factorial(n - k)) * 
					(*objects)[this->arg_indexes[3]].return_value(objects, args);
			}

			args->erase(args->begin() + var_index);
			return res;
		}
	}
	}
	if (this->type == ObjType::_Polynomial && args_results.size() > 3) {
		int var_index = args->size();
		args->push_back(Variable((*objects)[this->arg_indexes[0]].name));

		std::vector<std::vector<math::number>> matrix(args_results.size() - 3, std::vector<math::number>(args_results.size() - 2, 1));
		for (int i = 0; i < matrix.size(); i++) {
			for (int j = 1; j < matrix.size(); j++) {
				for (int k = j; k < matrix.size(); k++) {
					matrix[i][k] = matrix[i][k] * args_results[i + 3];
				}
			}
			(*args)[var_index].value = args_results[i + 3];
			matrix[i][matrix.size()] = (*objects)[this->arg_indexes[1]].return_value(objects, args);
		}

		for (int i = 0; i < matrix.size(); i++) {
			math::number multiple = matrix[i][i];
			for (int j = i; j < matrix[0].size(); j++)
				matrix[i][j] = matrix[i][j] / multiple;
			for (int j = 0; j < matrix.size(); j++) {
				if (i == j) continue;
				multiple = matrix[j][i];
				for (int k = 0; k < matrix[0].size(); k++)
					matrix[j][k] = matrix[j][k] - multiple * matrix[i][k];
			}
		}

		math::number
			var = 1,
			res = matrix[0][matrix.size()];
		for (int i = 1; i < matrix.size(); i++) {
			var = var * args_results[2];
			res = res + var * matrix[i][matrix.size()];
		}

		args->erase(args->begin() + var_index);
		return res;
	}
#if !infsimIsHere
	if (this->type == ObjType::_rand) {
		if (args_results.size() == 0) return math::rand(std::chrono::steady_clock::now().time_since_epoch().count() % 4096);
		return math::rand(0, args_results);
	}
#endif

	for (auto arg : *args)
		if (this->id == arg.id)
			return arg.value;
	for (auto fnc : functions_list)
		if (this->id == fnc.id) {
			for (int i = 0; i < fnc.args.size() && i < args_results.size(); i++)
				fnc.args[i].value = args_results[i];
			return fnc.return_value();
		}
	return 0;
}

math::number value(std::vector<Object> objects, std::vector<Variable> args) {
	return objects[0].return_value(&objects, &args);
}
math::number Function::return_value() {
	return value(this->objects, this->args);
}