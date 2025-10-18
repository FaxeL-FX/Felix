#pragma once
#include "../include.h"

#define infsimIsHere 0

namespace math {
	extern const long double pi, inf;

	bool sign(double);
	int E(double);
	unsigned long long M(double);
	double fract(double);

	//	complex
	extern struct complex_exponential;
	extern struct infsim;
	struct complex {
		long double R, i;
		complex() {
			this->R = 0;
			this->i = 0;
		}
		complex(long double R) {
			this->R = R;
			this->i = 0;
		}
		complex(long double R, long double i) {
			this->R = R;
			this->i = i;
		}
		complex(infsim x);
		bool isZero() {
			if (this->R == 0 && this->i == 0) return true;
			return false;
		}
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
	};
	extern const complex
		i,
		fctIntegralConstant;

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

	double factorial(double x);

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

	complex exp(complex x);
	complex ln(complex x);
	complex pow(complex x, complex y);
	complex sqrt(complex x);
	complex inv_sqrt(complex x);

	complex conjugate(complex x);

	complex cosh(complex x);
	complex sinh(complex x);
	complex coth(complex x);
	complex tanh(complex x);
	complex arccosh(complex x);
	complex arcsinh(complex x);
	complex arccoth(complex x);
	complex arctanh(complex x);

	complex cos(complex x);
	complex sin(complex x);
	complex cot(complex x);
	complex tan(complex x);
	complex arccos(complex x);
	complex arcsin(complex x);
	complex arccot(complex x);
	complex arctan(complex x);

	complex sin1(complex x);
	complex cos1(complex x);
	complex Binom(complex n, complex k);

	complex fct(complex x);
	complex fctIntegral(complex x, complex n);
	complex Harmonic(complex x);
	complex zeta(complex x);


	// infinitesimal (infsim)
	const unsigned int accuracy = 8u;
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
			for (higher; 0 <= higher && this->getNum(higher).isZero(); higher--) {}
			for (lower; lower < accuracy && this->getNum(lower).isZero(); lower++) {}
			for (int i = higher; i >= lower; i--) {
				std::string numstr = this->Pol[i].toString();
				if (numstr == "0") continue;
				str += "\n   + ";
				/**/ if (i == acch + 1)		str += "inf";
				else if (i == acch - 1)		str += "0";
				else if (i > acch)	str += "inf^" + std::to_string(i - acch);
				else if (i < acch)	str += "0^" + std::to_string(acch - i);
				else {
					str += "(" + numstr + ")";
					continue;
				}
				if (this->Pol[i].R > 0.999999 && this->Pol[i].R < 1.000001 &&
					this->Pol[i].i > -0.000001 && this->Pol[i].i < 0.000001) {
					continue;
				}
				str += "*(" + numstr + ")";
			}
			if (str == "") return "ConstZero";
			return str.substr(6);
		}
		std::string toStringSmall() {
			std::string str;
			int index = accuracy - 1;
			for (index; 0 <= index && this->getNum(index).isZero(); index--) {}
			if (index > acch)		return "inf";
			if (index < acch)		return "0";
			return this->Pol[index].toString();
		}
		std::string toStringLatex() {
			std::string str = "";
			int lower = 0, higher = accuracy - 1;
			for (higher; 0 <= higher && this->getNum(higher).isZero(); higher--) {}
			for (lower; lower < accuracy && this->getNum(lower).isZero(); lower++) {}
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
	extern const infsim infinity, zero, lnInf, fctIntegralConst;

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

	infsim cosh(infsim x);
	infsim sinh(infsim x);
	infsim coth(infsim x);
	infsim tanh(infsim x);
	infsim arccosh(infsim x);
	infsim arcsinh(infsim x);
	infsim arccoth(infsim x);
	infsim arctanh(infsim x);

	infsim cos(infsim x);
	infsim sin(infsim x);
	infsim cot(infsim x);
	infsim tan(infsim x);
	infsim arccos(infsim x);
	infsim arcsin(infsim x);
	infsim arccot(infsim x);
	infsim arctan(infsim x);

	infsim sin1(infsim x);
	infsim cos1(infsim x);
	infsim Binom(infsim n, infsim k);

	infsim fct(infsim x);
	infsim getFctIntegralConst();
	infsim fctIntegral(infsim x);
	infsim fctIntegral(infsim x, infsim y);
	infsim Harmonic(infsim x);
	infsim zeta(infsim x);

#if infsimIsHere
	typedef infsim number;
#else
	typedef complex number;
#endif
}