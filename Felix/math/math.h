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
	extern struct complex_exponential;
	struct fixed_point32 {
		int64_t num = 0;
		fixed_point32() {}
		fixed_point32(int n) {
			this->num = (long long)n << 32;
		}
		fixed_point32(long long n) {
			this->num = n << 32;
		}
		fixed_point32(long double n) {
			this->num = n * ((long long)1 << 32);
		}
		operator long double() {
			return (long double)this->num / ((long long)1 << 32);
		}
	};
	struct complex_linear {
		number R, i;
		complex_linear() {
			this->R = 0;
			this->i = 0;
		}
		complex_linear(number R) {
			this->R = R;
			this->i = 0;
		}
		complex_linear(number R, number i) {
			this->R = R;
			this->i = i;
		}
		complex_linear(complex_exponential x);
	};
	struct complex_exponential {
		long double r;
		fixed_point32 a;
		complex_exponential() {
			this->r = 0;
			this->a = 0;
		}
		complex_exponential(number x);
		complex_exponential(number r, number a) {
			this->r = r;
			this->a = a;
		}
		complex_exponential(complex_linear x);
	};
#if 1
	typedef complex_exponential complex;
#else
	typedef complex_linear complex;
#endif
	extern const complex_linear
		i,
		fctIntegralConstant;

	fixed_point32 operator+(fixed_point32, fixed_point32);
	fixed_point32 operator-(fixed_point32);
	fixed_point32 operator-(fixed_point32, fixed_point32);
	fixed_point32 operator*(fixed_point32, long double);
	fixed_point32 operator/(fixed_point32, long double);
	fixed_point32 operator%(fixed_point32, fixed_point32);
	fixed_point32 operator<<(fixed_point32, int);
	fixed_point32 operator>>(fixed_point32, int);
	bool operator==(fixed_point32, fixed_point32);

	complex_linear getFctIntegralConstant();

	complex_linear operator+(complex_linear, complex_linear);
	complex_linear operator-(complex_linear);
	complex_linear operator-(complex_linear, complex_linear);
	complex_linear operator*(complex_linear, complex_linear);
	complex_linear operator/(complex_linear, complex_linear);
	complex_linear operator%(complex_linear, complex_linear);
	complex_linear operator<<(complex_linear, int);
	complex_linear operator>>(complex_linear, int);
	bool operator==(complex_linear, complex_linear);

	complex_exponential operator+(complex_exponential, complex_exponential);
	complex_exponential operator-(complex_exponential);
	complex_exponential operator-(complex_exponential, complex_exponential);
	complex_exponential operator*(complex_exponential, complex_exponential);
	complex_exponential operator/(complex_exponential, complex_exponential);
	complex_exponential operator%(complex_exponential, complex_exponential);
	complex_exponential operator<<(complex_exponential, int);
	complex_exponential operator>>(complex_exponential, int);
	bool operator==(complex_exponential, complex_exponential);

	//	functions
	long double rand(int, std::vector<long double>);
	long double rand(int, std::vector<complex>);
	long double rand(int, complex_linear);
	long double rand(int);

	//	long double
	long double floor(long double);
	long double sign(long double);

	long double exp(long double);
	long double ln(long double);

	long double sqrt(long double);
	long double inv_sqrt(long double);

	long double cos(long double);
	long double sin(long double);
	long double arccos(long double);

	//	complex_linear
	std::string toString(complex_linear);

	complex_linear floor(complex_linear);
	long double abs(complex_linear);
	long double inv_abs(complex_linear);
	complex_linear normalize(complex_linear);
	complex_linear mul_i(complex_linear);
	long double arg(complex_linear);
	bool exist(complex_linear);

	complex_linear exp(complex_linear);
	complex_linear ln(complex_linear);
	complex_linear pow(complex_linear, complex_linear);
	complex_linear sqrt(complex_linear);
	complex_linear inv_sqrt(complex_linear);

	//	complex_exponential
	std::string toString(complex_exponential);
	complex_exponential floor(complex_exponential);
	long double abs(complex_exponential);
	long double inv_abs(complex_exponential);
	complex_exponential normalize(complex_exponential);
	complex_exponential mul_i(complex_exponential);
	long double arg(complex_exponential);
	bool exist(complex_exponential);

	complex_exponential exp(complex_exponential);
	complex_exponential ln(complex_exponential);
	complex_exponential pow(complex_exponential, complex_exponential);
	complex_exponential sqrt(complex_exponential);
	complex_exponential inv_sqrt(complex_exponential);

	//	complex
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
	complex fctIntegral(complex, complex);
}