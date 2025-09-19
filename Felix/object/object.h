#pragma once
#include "../include.h"
#include "../math/math.h"

enum ObjType {
	_Error, _Default, _Const, _rand,

	_Sum_List, _Product_List,
	_Sum_FOR, _Product_FOR,
	_flip_sign, _inverse,

	_exp, _ln,
	_Harmonic,

	_Integral_Defined, _Integral_Undefined,
	_Return,



	_add, _dif, _mul, _div, _pow, _fct, _mod, _uMinus,

	_sqrt, _inv_sqrt,
	_root, _log,

	_cos, _arccos,
	_sin, _arcsin,
	_cot, _arccot,
	_tan, _arctan,
	_cosh, _arccosh,
	_sinh, _arcsinh,
	_coth, _arccoth,
	_tanh, _arctanh,
	_sin1, _cos1,

	_Sum, _ForwardDifference,
	_Product,
	_Integral, _Derivative, _IntegralAlongExp,
	_Polynomial,

	_abs, _inv_abs, _arg, _sign, _Re, _Im, _floor, _ceil, _round,

	_exist, _grid,

	_gamma, _fctIntegral, _zeta,
};

struct Variable {
	int id = -1;
	std::string name;
	math::number value = 0;

	Variable(std::string name) {
		this->name = name;
		this->id = 0;
		for (auto c : this->name)
			this->id = (this->id << 1) + (int)c;
	}
	Variable(std::string name, math::number value) {
		this->name = name;
		this->value = value;
		this->id = 0;
		for (auto c : this->name)
			this->id = (this->id << 1) + (int)c;
	}
	Variable(int id, std::string name, math::number value) {
		this->id = id;
		this->name = name;
		this->value = value;
	}
	Variable() {}
};
struct Object {
	int id;
	std::string name;
	ObjType type;
	std::vector<int> arg_indexes;
	math::number value;

	math::number return_value(std::vector<Object>* objects, std::vector<Variable>* args);

	Object(std::string name, ObjType type, std::vector<int> arg_indexes, math::number value) {
		this->name = name;
		this->type = type;
		this->arg_indexes = arg_indexes;
		this->value = value;
		this->id = 0;
		for (auto c : this->name)
			this->id = (this->id << 1) + (int)c;
	}
	Object(std::string name, std::vector<int> arg_indexes, math::number value) {
		this->name = name;
		this->arg_indexes = arg_indexes;
		this->value = value;
		this->id = 0;
		for (auto c : this->name)
			this->id = (this->id << 1) + (int)c;
		this->type = ObjType::_Default;
	}
	Object(ObjType type, std::vector<int> arg_indexes) {
		this->type = type;
		this->arg_indexes = arg_indexes;
	}
	Object(std::string name) {
		this->name = name;
		this->id = 0;
		for (auto c : this->name)
			this->id = (this->id << 1) + (int)c;
		this->type = ObjType::_Default;
	}
	Object() {
		this->id = -1;
	}
};
struct Function {
	int id = -1;
	std::string name;
	std::vector<Variable> args;
	std::vector<Object> objects;

	math::number return_value();

	Function() {}
	Function(int id, std::string name, std::vector<Variable> args, std::vector<Object> objects) {
		this->id = id;
		this->name = name;
		this->args = args;
		this->objects = objects;
	}
};

extern std::vector<Function> functions_list;

ObjType nameToType(std::string);

std::vector<Object> parse_expr(std::string);
std::string parse_token(std::string);
void parse_obj(std::vector<Object> *objects, std::vector<std::string> tokens, int begin, int end);

int brackets(std::vector<std::string>, int, bool);
char antibr(char);
int priority(std::string);
int find_token(std::vector<std::string>, int, bool);
bool need_mul(std::string, std::string);

math::number value(std::vector<Object>, std::vector<Variable>);