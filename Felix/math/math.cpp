#include "math.h"
namespace math {
	//	number
	const number pi = 3.1415926535897932384626433832795028841971693993751058209749445923;

	//number operator+(number x, number y) { return x + y; }
	//number operator-(number x) { return -x; }
	//number operator-(number x, number y) { return x - y; }
	//number operator*(number x, number y) { return x * y; }
	//number operator/(number x, number y) { return x / y; }
	//number operator%(number x, number y) { return x - y * floor(x / y); }
	//number operator<<(number x, int n) { return 0; }
	//number operator>>(number x, int n) { return x << -n; }

	bool sign(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return i >> 63;
	}
	int E(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return (int)((i >> 52) % 2048) - 1023;
	}
	unsigned long long M(double x) {
		unsigned long long i = *(unsigned long long*) & x;
		return i % ((unsigned long long)1 << 52);
	}
	double fract(double x) {
		return (double)M(x) / ((unsigned long long)1 << 52);
	}

	//	complex
	const complex i(0, 1);

	complex operator+(complex x, complex y) { return complex(x.R + y.R, x.i + y.i); }
	complex operator-(complex x) { return complex(-x.R, -x.i); }
	complex operator-(complex x, complex y) { return complex(x.R - y.R, x.i - y.i); }
	complex operator*(complex x, complex y) { return complex(x.R * y.R - x.i * y.i, x.R * y.i + x.i * y.R); }
	complex operator/(complex x, complex y) {
		number denominator = 1 / (y.R * y.R + y.i * y.i);
		return complex((x.R * y.R + x.i * y.i) * denominator, (y.R * x.i - x.R * y.i) * denominator);
	}
	complex operator%(complex x, complex y) { return x - y * floor(x / y); }
	complex operator<<(complex x, int n) {
		if (n < 0) {
			for (int i = 0; i < -n / 64; i++) x = x / ((unsigned long long)1 << 63) / 2;
			return x / ((unsigned long long)1 << -n % 64);
		}
		for (int i = 0; i < n / 64; i++) x = x * ((unsigned long long)1 << 63) * 2;
		return x * ((unsigned long long)1 << n % 64);
	}																						//	need to change for "flex_float" instead "long double"
	complex operator>>(complex x, int n) { return x << -n; }

	//	functions
	number rand(int seed, std::vector<number> v) {
		number res = cos(seed);
		for (int i = 1; i <= 16; i++)
			for (auto n : v)
				res = res + sin(n) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}																						//	need to change for "flex_float" instead "long double"
	number rand(int seed, complex z) {
		std::vector<number> v = { z.R, z.i };
		return rand(seed, v);
	}
	number rand(int n) {
		std::vector<number> v = { (number)n };
		return rand(n, v);
	}

	//	for number
	number floor(number x) {
		return (long long)x;
	}																						//	need to change for "flex_float" instead "long double"
	number sign(number x) {
		if (x < 0) return -1;
		return 1;
	}

	number exp(number x) {
		number res = x / ((unsigned long long)1 << 32);
		if (res < 0) res = -res;
		for (int i = 0; i < 32; i++) res = 2 * res + res * res;
		if (x < 0) return 1 / (res + 1);
		return res + 1;
	}																						//	need to change for "flex_float" instead "long double"
	number log(number x) {
		number res = 0, part = 1;
		for (int k = 1; k < 32; k++) {
			part = part * (1 - x);
			res = res + ((1 - part) / k) / ((unsigned long long)1 << k);
		}
		return res;
	}																						//	need to change for "flex_float" instead "long double"
	number ln(number x) { return log(fract(x)) + 0.6931471805599453 * E(x); }

	number sqrt(number x) { return exp(0.5 * ln(x)); }
	number inv_sqrt(number x) { return exp(-0.5 * ln(x)); }

	number cos(number x) {
		number res = x / ((unsigned long long)1 << 16);
		res = 2 - res * res;
		for (int i = 0; i < 16; i++) res = res * res - 2;
		return res * 0.5;
	}																						//	need to change for "flex_float" instead "long double"
	number sin(number x) { return cos(x - 1.57079632679); }
	number arccos(number x) {
		if (x == 1) return 0;
		if (x < 0) return pi - arccos(-x);
		number res = sqrt((x + 1) * 2);
		for (int i = 1; i < 12; i++) res = sqrt(2 + res);
		return sqrt(2 - res) * ((unsigned long long)1 << 12);
	}																						//	need to change for "flex_float" instead "long double"

	//	for complex
	complex floor(complex x) { return complex(floor(x.R), floor(x.i)); }
	number abs(complex x) {
		if (x.i == 0) return sign(x.R) * x.R;
		if (x.R == 0) return sign(x.i) * x.i;
		return sqrt(x.R * x.R + x.i * x.i);
	}
	number inv_abs(complex x) { return inv_sqrt(x.R * x.R + x.i * x.i); }
	complex normalize(complex x) {
		if (x.i == 0)
			if (x.R < 0)	return -1;
			else			return 1;
		if (x.R == 0)
			if (x.i < 0)	return -i;
			else			return i;
		return x * inv_sqrt(x.R * x.R + x.i * x.i);
	}
	complex mul_i(complex x) { return complex(-x.i, x.R); }

	complex exp(complex x) {
		complex res = x >> 64;
		if (res.R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.R < 0) return 1 / (res + 1);
		return res + 1;
	}																						//	need to change for "flex_float" instead "long double"
	complex ln(complex x) {
		number cosine;
		if (x.i == 0)
			if (x.R < 0)	cosine = -1;
			else			cosine = 1;
		else				cosine = x.R * inv_abs(x);
		if (x.i < 0) return complex(ln(abs(x)), -arccos(cosine));
		return complex(ln(abs(x)), arccos(cosine));
	}																						//	need to change for "flex_float" instead "long double"

	complex pow(complex x, complex y) { return exp(y * ln(x)); }
	complex sqrt(complex x) { return exp(0.5 * ln(x)); }
	complex inv_sqrt(complex x) { return exp(-0.5 * ln(x)); }

	complex cosh(complex x) { return (exp(x) + exp(-x)) * 0.5; }
	complex sinh(complex x) { return (exp(x) - exp(-x)) * complex(0, 0.5); }
	complex coth(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	complex tanh(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	complex arccosh(complex x) { return ln(x + sqrt(x * x - 1)); }
	complex arcsinh(complex x) { return ln(x + sqrt(x * x + 1)); }
	complex arccoth(complex x) { return arctanh(1 / x); }
	complex arctanh(complex x) { return -ln((1 - x) * inv_sqrt(1 - x * x)); }

	complex cos(complex x) { return        cosh(mul_i(x)); }
	complex sin(complex x) { return -mul_i(sinh(mul_i(x))); }
	complex cot(complex x) { return  mul_i(coth(mul_i(x))); }
	complex tan(complex x) { return -mul_i(tanh(mul_i(x))); }
	complex arccos(complex x) { return -mul_i(arccosh(x)); }
	complex arcsin(complex x) { return -mul_i(arcsinh(mul_i(x))); }
	complex arccot(complex x) { return  mul_i(arccoth(mul_i(x))); }
	complex arctan(complex x) { return -mul_i(arctanh(mul_i(x))); }

	complex fct(complex x) {
		if (x.R < -0.5) return -pi / (sin(pi * x) * fct(-1 - x));
		const long double n = 32.5;
		complex res = exp(-x + (x + n) * ln(x + n) - n * ln(n));
		for (int i = 1; i < n; i++) res = res * i / (x + i);
		return res;
	}
}