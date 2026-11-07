#pragma once
#include "../include.h"
#include "fixed_point.h"

namespace math {
	//using real = math::fixed_point<32>;
	using real = long double;
	extern const real pi, inf;

	inline std::string to_string(long double x, size_t digits) {
		return std::format("{:.{}f}", x, digits);
	}
	inline std::string to_string(long long x, size_t digits) {
		return std::format("{}", x);
	}
	template<size_t bit>
	inline std::string to_string(math::fixed_point<bit> x, size_t digits) {
		return x.toString(digits);
	}

	extern struct polar;

	//	complex
	struct complex {
		real R, i;

		complex() {
			this->R = 0;
			this->i = 0;
		}
		complex(real R) {
			this->R = R;
			this->i = 0;
		}
		complex(real R, real i) {
			this->R = R;
			this->i = i;
		}
		complex(polar x);

		std::string toString(size_t digits = 20) {
			real R = this->R, i = this->i;

			if (i == (real)0) return to_string(this->R, digits);
			if (R == (real)0) return to_string(this->i, digits) + "i";

			if (i < (real)0) return to_string(this->R, digits) + to_string(this->i, digits) + "i";
			return to_string(this->R, digits) + "+" + to_string(this->i, digits) + "i";
		}
	};
	extern const complex i;

	complex operator+(complex, complex);
	complex operator-(complex);
	complex operator-(complex, complex);
	complex operator*(complex, complex);
	complex operator/(complex, complex);
	complex operator%(complex, complex);
	bool operator==(complex, complex);
	bool operator!=(complex, complex);


	struct polar {
		real radius, angle;

		polar() {
			this->radius = 0;
			this->angle = 0;
		}
		polar(real x) {
			if (x < 0) {
				this->radius = -x;
				this->angle = pi;
			}
			else {
				this->radius = x;
				this->angle = 0;
			}
		}
		polar(real radius, real angle) {
			this->radius = radius;
			this->angle = angle;
		}
		polar(complex x);

		std::string toString(size_t digits = 20) {
			return complex(*this).toString(digits);
		}
	};

	polar operator-(polar);
	polar operator*(polar, polar);
	polar operator/(polar, polar);
	bool operator==(polar, polar);
	bool operator!=(polar, polar);



	//	FUNCTIONS

	//	real
	real rand(int, std::vector<complex>);
	real rand(int);

	real factorial(real x);
	real Binom(real n, real k);
	extern std::vector<real> zbfValues;
	real zbf(int x);
	extern std::vector<real> zetaZeros;
	real zetaZero(int n);

	//	complex
	complex floor(complex x);
	complex ceil(complex x);
	real abs(complex x);
	complex normalize(complex x);
	complex mul_i(complex x);
	real arg(complex x);
	bool exist(complex x);
	real Re(complex x);
	real Im(complex x);
	real grid(complex x);
	real axis(complex x);

	complex exp(complex x);
	complex ln(complex x);
	complex pow(complex x, complex y);
	complex sqrt(complex x);
	complex inv_sqrt(complex x);

	inline complex cosh(complex x) { return (exp(x) + exp(-x)) * (real)0.5; }
	inline complex sinh(complex x) { return (exp(x) - exp(-x)) * (real)0.5; }
	inline complex coth(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	inline complex tanh(complex x) {
		complex positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	inline complex arccosh(complex x) { return ln(x + sqrt(x * x - (real)1)); }
	inline complex arcsinh(complex x) { return ln(x + sqrt(x * x + (real)1)); }
	inline complex arccoth(complex x) { return (real)-0.5 * ln((x - (real)1) / (x + (real)1)); }
	inline complex arctanh(complex x) { return (real)-0.5 * ln(((real)1 - x) / ((real)1 + x)); }

	inline complex cos(complex x) { return        cosh(mul_i(x)); }
	inline complex sin(complex x) { return -mul_i(sinh(mul_i(x))); }
	inline complex cot(complex x) { return  mul_i(coth(mul_i(x))); }
	inline complex tan(complex x) { return -mul_i(tanh(mul_i(x))); }
	inline complex arccos(complex x) { return -mul_i(arccosh(x)); }
	inline complex arcsin(complex x) { return -mul_i(arcsinh(mul_i(x))); }
	inline complex arccot(complex x) { return  mul_i(arccoth(mul_i(x))); }
	inline complex arctan(complex x) { return -mul_i(arctanh(mul_i(x))); }

	complex sin1(complex x);
	complex cos1(complex x);
	complex Binom(complex n, complex k);

	complex fct(complex x);
	complex gamma(complex x);
	complex fctIntegral(complex x, complex n);
	complex Harmonic(complex x);
	complex zeta(complex x);
	complex zetaByFct(complex x);

	complex USumN(complex x, complex n);

	inline complex UIntN(complex x, complex n) {
		if (n == complex{ -1 }) return ln(x);
		return pow(x, n + complex{ 1 }) / (n + complex{ 1 });
	}
	inline complex Ei(complex x) {
		complex res{ 0.5772156649 }, u{ 1 };
		for (int k = 1; k < 256; k++) {
			u = u * x / (real)k;
			res = res + u / (real)k;
		}
		return res;
	}
	inline complex Li(complex x) { return ln(ln(x)) + Ei(ln(x)); }

	complex Integral(complex(*f)(complex), complex a, complex b);
	complex Integral(complex(*f)(complex), complex(*path)(real));
	complex Derivative(complex(*f)(complex), complex x);

	// polar
	polar pow(polar x, polar y);
	polar exp(polar x);
	polar ln(polar x);
	polar sqrt(polar x);



	typedef complex number;
}