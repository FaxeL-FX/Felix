#pragma once
#include "../include.h"

namespace math {
	//	number
	//	flexible integer
	struct flex_uint {
		std::vector<unsigned long long> flex;

		flex_uint() {}
		flex_uint(unsigned long long) {}
		flex_uint(std::vector<unsigned long long> flex) {
			this->flex = flex;
		}
		flex_uint(std::string) {}
		std::string toString() {
			return "";
		}
	};
	//	signed
	struct flex_int {
		std::vector<unsigned long long> flex;
		bool sign = false;

		flex_int() {}
		flex_int(long long) {}
		flex_int(std::string) {}
		std::string toString() {
			std::string s = "";
			if (sign) s = "-";
			return s + flex_uint(this->flex).toString();
		}
	};
	//	unsigned fraction
	struct flex_ufract {
		std::vector<unsigned long long> flex;

		flex_ufract() {}
		flex_ufract(unsigned long long) {}
		flex_ufract(std::string) {}
		std::string toString() {
			return "";
		}
	};
	//	flexible float point
	struct flex_float {
		bool		sign = false;
		flex_int	exponent = 0;
		flex_ufract	mantissa = 0;

		flex_float() {}
		flex_float(double) {}
		flex_float(std::string) {}
		std::string toString() {
			return "";
		}
	};

	//	def number
	typedef long double number;
	extern const number pi;

	bool sign(double);
	int E(double);
	unsigned long long M(double);
	double fract(double);

	//	complex
	struct complex {
		number R, i;
		complex() { R = 0; i = 0; }
		complex(number R) { this->R = R; i = 0; }
		complex(number R, number i) { this->R = R; this->i = i; }
	};
	extern const complex i, fctIntegralConstant;

	complex getFctIntegralConstant();

	complex operator+(complex, complex);
	complex operator-(complex);
	complex operator-(complex, complex);
	complex operator*(complex, complex);
	complex operator/(complex, complex);
	complex operator%(complex, complex);
	complex operator<<(complex, int);
	complex operator>>(complex, int);

	//	functions
	long double rand(int, std::vector<long double>);
	long double rand(int, std::vector<complex>);
	long double rand(int, complex);
	long double rand(int);

	//	for number
	long double floor(long double);
	long double sign(long double);

	long double exp(long double);
	long double ln(long double);

	long double sqrt(long double);
	long double inv_sqrt(long double);

	long double cos(long double);
	long double sin(long double);
	long double arccos(long double);

	//	for complex
	bool operator==(complex, complex);
	std::string toString(complex);

	complex floor(complex);
	long double abs(complex);
	long double inv_abs(complex);
	complex normalize(complex);
	complex mul_i(complex);
	long double arg(complex);

	complex exp(complex);
	complex ln(complex);

	complex pow(complex, complex);
	complex sqrt(complex);
	complex inv_sqrt(complex);

	complex cosh(complex);
	complex sinh(complex);
	complex coth(complex);
	complex tanh(complex);
	complex arccosh(complex);
	complex arcsinh(complex);
	complex arccoth(complex);
	complex arctanh(complex);

	complex cos(complex);
	complex sin(complex);
	complex cot(complex);
	complex tan(complex);
	complex arccos(complex);
	complex arcsin(complex);
	complex arccot(complex);
	complex arctan(complex);

	complex fct(complex);
	complex fctIntegral(complex);
}