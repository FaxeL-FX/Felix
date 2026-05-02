#pragma once
#include "../include.h"

#ifndef infsimIsHere
#define infsimIsHere 0
#endif

namespace math {
	extern const double pi, inf;

	bool sign(double);
	int E(double);
	unsigned long long M(double);
	double fract(double);

	//	complex
	extern struct complex_exponential;
	extern struct infsim;
	using real = long double;
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
		complex(infsim x);

		std::string toString() {
			double R = this->R, i = this->i;
			if (-0.000001 < this->R && this->R < 0.000001) R = 0.0;
			if (-0.000001 < this->i && this->i < 0.000001) i = 0.0;

			if (R == 0.0 && this->i == 0.0) return "0";
			if (R == 0.0) return std::to_string(this->i) + "i";
			if (i == 0.0) return std::to_string(this->R);

			if (i < 0) return std::to_string(this->R) + std::to_string(this->i) + "i";
			return std::to_string(this->R) + "+" + std::to_string(this->i) + "i";
		}

		void operator +=(complex& x) {
			this->R += x.R;
			this->i += x.i;
		}
		void operator -=(complex& x) {
			this->R -= x.R;
			this->i -= x.i;
		}
		void operator *=(complex& y) {
			real 
				new_R = this->R * y.R - this->i * y.i,
				new_i = this->R * y.i + this->i * y.R;
			this->R = new_R;
			this->i = new_i;
		}
		void inverse() {
			real denominator = this->R * this->R + this->i * this->i;
			this->R /= denominator;
			this->i /= denominator;
			this->i = -this->i;
		}
		void operator /=(complex y) {
			y.inverse();
			*this *= y;
		}
	};
	extern const complex i;

	complex getFctIntegralConstant();

	complex operator+(complex, complex);
	complex operator-(complex);
	complex operator-(complex, complex);
	complex operator*(complex, complex);
	complex operator/(complex, complex);
	complex operator%(complex, complex);
	complex operator<<(complex, int);
	complex operator>>(complex, int);
	bool operator==(complex, complex);
	bool operator!=(complex, complex);

	//	functions
	long double rand(int, std::vector<long double>);
	long double rand(int, std::vector<complex>);
	long double rand(int, complex);
	long double rand(int);

	//	long double
	long double floor(long double);
	long double ceil(long double);
	long double sign(long double);
	long double round(long double);

	long double exp(long double);
	long double ln(long double);
	long double pow(long double, int);

	long double sqrt(long double);
	long double inv_sqrt(long double);

	long double cos(long double);
	long double sin(long double);
	long double arccos(long double);

	long double factorial(long double x);
	long double Binom(long double n, long double k);
	extern std::vector<long double> zbfValues;
	long double zbf(int x);
	extern std::vector<long double> zetaZeros;
	long double zetaZero(int n);

	//	complex
	complex floor(complex);
	complex ceil(complex);
	long double abs(complex);
	long double inv_abs(complex);
	complex normalize(complex);
	complex mul_i(complex);
	long double arg(complex);
	bool exist(complex);
	complex Re(complex x);
	complex Im(complex x);
	complex grid(complex x);
	complex axis(complex x);

	complex exp(complex x);
	complex ln(complex x);
	complex pow(complex x, complex y);
	complex sqrt(complex x);
	complex inv_sqrt(complex x);

	complex conjugate(complex x);

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


	// infinitesimal (infsim)
	const unsigned int accuracy = 32u;
	const int acch = accuracy >> 1;

	struct infsim {
		std::vector<complex> Pol;
		infsim() { this->Pol.resize(accuracy); }
		infsim(complex x) {
			this->Pol.resize(accuracy);
			Pol[acch] = x;
		}
		infsim(double x) : infsim(complex(x)) {}
		infsim(int index, complex x) {
			this->Pol.resize(accuracy);
			if (index < 0) return;
			if (index >= accuracy) return;
			Pol[index] = x;
		}
		complex getNum(int index) {
			if (index < 0) return 0;
			if (index >= accuracy) return 0;
			return Pol[index];
		}
		void setNum(int index, complex x) {
			if (index < 0) return;
			if (index >= accuracy) return;
			Pol[index] = x;
		}
		std::string toString() {
			std::string str = "";
			int lower = 0, higher = accuracy - 1;
			for (; 0 <= higher && this->getNum(higher).isZero(); higher--) {}
			for (; lower < accuracy && this->getNum(lower).isZero(); lower++) {}
			for (int i = higher; i >= lower; i--) {
				std::string numstr = this->Pol[i].toString();
				if (numstr == "0") continue;
				str += "\n  + (" + numstr + ")";
				/**/ if (i == acch + 1)		str += " * inf";
				else if (i == acch - 1)		str += " * 0";
				else if (i > acch)	str += " * inf^" + std::to_string(i - acch);
				else if (i < acch)	str += " * 0^" + std::to_string(acch - i);
			}
			if (str == "") return "ConstZero";
			return str.substr(5);
		}
		std::string toStringSmall() {
			std::string str;
			int index = accuracy - 1;
			for (; 0 <= index && this->getNum(index).isZero(); index--) {}
			if (index > acch)		return "inf";
			if (index < acch)		return "0";
			return this->Pol[index].toString();
		}
		std::string toStringLatex() {
			std::string str = "";
			int lower = 0, higher = accuracy - 1;
			for (; 0 <= higher && this->getNum(higher).isZero(); higher--) {}
			for (; lower < accuracy && this->getNum(lower).isZero(); lower++) {}
			for (int i = higher; i >= lower; i--) {
				if (i < higher)		str += "+";
				/**/ if (i == acch)		str += "1\\cdot";
				else if (i == acch + 1)		str += "\\infty\\cdot";
				else if (i == acch - 1)		str += "0\\cdot";
				else if (i > acch)	str += "\\infty^" + std::to_string(i - acch) + "\\cdot";
				else				str += "0^" + std::to_string(acch - i) + "\\cdot";
				str += "(" + this->Pol[i].toString() + ")";
			}
			if (str == "") str = "0_C";
			return str;
		}
	};
	extern const infsim infinity, zero, lnInf;

	infsim operator+(infsim x, infsim y);
	infsim operator-(infsim x);
	infsim operator-(infsim x, infsim y);
	infsim operator*(infsim x, infsim y);
	infsim inv(infsim x);
	infsim operator/(infsim x, infsim y);
	infsim mul(infsim x, complex y);
	infsim div(infsim x, complex y);
	infsim operator<<(infsim x, int n);
	infsim operator>>(infsim x, int n);


	infsim Re(infsim x);
	infsim Im(infsim x);
	complex grid(infsim x);
	infsim mul_i(infsim x);
	infsim floor(infsim x);
	infsim ceil(infsim x);

	infsim exp(infsim x);
	infsim ln_sum(infsim x);
	infsim ln_integral(infsim x);
	infsim ln(infsim x);
	infsim pow(infsim x, infsim y);
	infsim sqrt(infsim x);
	infsim inv_sqrt(infsim x);

	infsim sin1(infsim x);
	infsim cos1(infsim x);
	infsim Binom(infsim n, infsim k);

	infsim fct(infsim x);
	infsim inv_fct(infsim x);
	infsim gamma(infsim x);
	infsim fctIntegral(infsim x, infsim y);
	infsim Harmonic(infsim x);
	infsim zeta(infsim x);

	infsim USumN(infsim x, infsim n);

#if infsimIsHere
	typedef infsim number;
#else
	typedef complex number;
#endif

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