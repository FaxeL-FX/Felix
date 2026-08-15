#pragma once
#include "../include.h"

namespace math {
	using real = long double;
	extern const real pi, inf;

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

		std::string toString() {
			double R = this->R, i = this->i;
			
			if (i == 0.0) return std::format("{:.20f}", this->R);
			if (R == 0.0) return std::format("{:.20f}", this->i) + "i";

			if (i < 0) return std::format("{:.20f}", this->R) + std::format("{:.20f}", this->i) + "i";
			return std::format("{:.20f}", this->R) + "+" + std::format("{:.20f}", this->i) + "i";
		}

		complex& operator +=(complex& x) {
			this->R += x.R;
			this->i += x.i;
			return *this;
		}
		complex& operator -=(complex& x) {
			this->R -= x.R;
			this->i -= x.i;
			return *this;
		}
		complex& operator *=(complex& y) {
			real 
				new_R = this->R * y.R - this->i * y.i,
				new_i = this->R * y.i + this->i * y.R;
			this->R = new_R;
			this->i = new_i;
			return *this;
		}
		void inverse() {
			real denominator = this->R * this->R + this->i * this->i;
			this->R /= denominator;
			this->i /= denominator;
			this->i = -this->i;
		}
		complex& operator /=(complex y) {
			y.inverse();
			*this *= y;
			return *this;
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

	complex sin1(complex x);
	complex cos1(complex x);
	complex Binom(complex n, complex k);

	complex fct(complex x);
	complex inv_fct(complex x);
	complex gamma(complex x);
	complex fctIntegral(complex x, complex n);
	complex Harmonic(complex x);
	complex zeta(complex x);
	complex zetaByFct(complex x);

	complex USumN(complex x, complex n);

	complex Integral(complex(*f)(complex), complex a, complex b);
	complex Integral(complex(*f)(complex), complex(*path)(real));
	complex Derivative(complex(*f)(complex), complex x);


	typedef complex number;


	template <typename NumT> NumT cosh(NumT x) { return (exp(x) + exp(-x)) * 0.5; }
	template <typename NumT> NumT sinh(NumT x) { return (exp(x) - exp(-x)) * 0.5; }
	template <typename NumT> NumT coth(NumT x) {
		NumT positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	template <typename NumT> NumT tanh(NumT x) {
		NumT positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	template <typename NumT> NumT arccosh(NumT x) { return ln(x + sqrt(x * x - 1)); }
	template <typename NumT> NumT arcsinh(NumT x) { return ln(x + sqrt(x * x + 1)); }
	template <typename NumT> NumT arccoth(NumT x) { return -0.5 * ln((x - 1) / (x + 1)); }
	template <typename NumT> NumT arctanh(NumT x) { return -0.5 * ln((1 - x) / (1 + x)); }

	template <typename NumT> NumT cos(NumT x) { return        cosh(mul_i(x)); }
	template <typename NumT> NumT sin(NumT x) { return -mul_i(sinh(mul_i(x))); }
	template <typename NumT> NumT cot(NumT x) { return  mul_i(coth(mul_i(x))); }
	template <typename NumT> NumT tan(NumT x) { return -mul_i(tanh(mul_i(x))); }
	template <typename NumT> NumT arccos(NumT x) { return -mul_i(arccosh(x)); }
	template <typename NumT> NumT arcsin(NumT x) { return -mul_i(arcsinh(mul_i(x))); }
	template <typename NumT> NumT arccot(NumT x) { return  mul_i(arccoth(mul_i(x))); }
	template <typename NumT> NumT arctan(NumT x) { return -mul_i(arctanh(mul_i(x))); }

	template <typename NumT> NumT UIntN(NumT x, NumT n) {
		if (n == NumT{ -1 }) return ln(x);
		return pow(x, n + 1) / (n + 1);
	}
	template <typename NumT> NumT Ei(NumT x) {
		NumT res = 0.5772156649, u = 1;
		for (int k = 1; k < 256; k++) {
			u = u * x / k;
			res = res + u / k;
		}
		return res;
	}
	template <typename NumT> NumT Li(NumT x) { return ln(ln(x)) + Ei(ln(x)); }
}