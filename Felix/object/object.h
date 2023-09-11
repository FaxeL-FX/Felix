#pragma once
#include "../include.h"
#include "../math/math.h"

namespace mobj {
	enum mathObjs {
		_Error, _Default, _Const, _Rand,

		_add, _dif, _mul, _div, _pow, _fct, _mod, _uMinus,

		abs, inv_abs, arg, Re, Im, exist, floor, round, grid,

		exp, ln,
		sqrt, inv_sqrt,
		root, log,

		cos, arccos,
		sin, arcsin,
		cot, arccot,
		tan, arctan,
		cosh, arccosh,
		sinh, arcsinh,
		coth, arccoth,
		tanh, arctanh,

		gamma,

		Sum, Product, Return,
		Integral, Derivative, IntegralAlongExp,
	};
	struct mathObject {
		mathObjs type;
		std::vector<std::string> names;
		math::complex(*function)(math::complex, ...);
		int args_number;
	};
}

struct Variable {
	std::string name;
	math::complex value;

	Variable(std::string name, math::complex value) {
		this->name = name;
		this->value = value;
	}
};
struct Object {
	std::string name;
	mobj::mathObjs type;
	std::vector<int> arg_indexes;
	math::complex value;

	math::complex return_value(std::vector<Object>* objects, std::vector<Variable>* args);

	Object(std::string name, std::vector<int> arg_indexes, math::complex value) {
		this->name = name;
		this->arg_indexes = arg_indexes;
		this->value = value;
	}
	Object(std::string name) {
		this->name = name;
	}
	Object() {}
};
struct Function {
	std::string name;
	std::vector<Variable> args;
	std::vector<Object> objects;

	math::complex return_value() {
		if (objects.size() < 1) return 0;
		std::vector<Object>* obj = &objects;
		std::vector<Variable>* ar = &args;
		return objects[0].return_value(obj, ar);
	}

	Function() {}
	Function(std::string name, std::vector<Variable> args, std::vector<Object> objects) {
		this->name = name;
		this->args = args;
		this->objects = objects;
	}
};

extern std::vector<Function> functions_list;

std::string toString(std::vector<Object>, int);
mobj::mathObjs nameToType(std::string);

std::vector<Object> parse_expr(std::string);
std::string parse_token(std::string);
void parse_obj(std::vector<Object>*, std::vector<std::string>, int, int);

int brackets(std::vector<std::string>, int, bool);
char antibr(char);
int priority(std::string);
int find_token(std::vector<std::string>, int, bool);
bool need_mul(std::string, std::string);

math::complex value(std::vector<Object>, std::vector<Variable>);