#pragma once
#include "../include.h"
#include "../math/math.h"

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
		std::vector<Object>* obj = &objects;
		std::vector<Variable>* ar = &args;
		return objects[0].return_value(obj, ar);
	}

	Function(std::string name, std::vector<Variable> args, std::vector<Object> objects) {
		this->name = name;
		this->args = args;
		this->objects = objects;
	}
};

std::vector<Object> parse_expr(std::string);
std::string parse_token(std::string);
void parse_obj(std::vector<Object>*, std::vector<std::string>, int, int);

int brackets(std::vector<std::string>, int, bool);
char antibr(char);
int priority(std::string);
int find_token(std::vector<std::string>, int, bool);
bool need_mul(std::string, std::string);