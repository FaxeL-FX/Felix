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

	extern const long double pi, inf;

	bool sign(double);
	int E(double);
	unsigned long long M(double);
	double fract(double);

	//	complex
	extern struct complex_exponential;
	extern struct infsim;
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
		long double R, i;
		complex_linear() {
			this->R = 0;
			this->i = 0;
		}
		complex_linear(long double R) {
			this->R = R;
			this->i = 0;
		}
		complex_linear(long double R, long double i) {
			this->R = R;
			this->i = i;
		}
		complex_linear(complex_exponential x);
		complex_linear(infsim x);
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
	struct complex_exponential {
		long double r, a;
		complex_exponential() {
			this->r = 0;
			this->a = 0;
		}
		complex_exponential(long double x);
		complex_exponential(long double r, long double a) {
			this->r = r;
			this->a = a;
		}
		complex_exponential(complex_linear x);
	};
#if 0
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
	fixed_point32 operator*(fixed_point32, fixed_point32);
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
	long double ceil(long double);
	long double sign(long double);
	long double round(long double);

	long double exp(long double);
	long double ln(long double);

	long double sqrt(long double);
	long double inv_sqrt(long double);

	long double cos(long double);
	long double sin(long double);
	long double arccos(long double);

	double factorial(double x);

	//	complex_linear
	std::string toString(complex_linear);

	complex_linear floor(complex_linear);
	complex_linear ceil(complex_linear);
	long double abs(complex_linear);
	long double inv_abs(complex_linear);
	complex_linear normalize(complex_linear);
	complex_linear mul_i(complex_linear);
	long double arg(complex_linear);
	bool exist(complex_linear);
	complex_linear Re(complex_linear x);
	complex_linear Im(complex_linear x);
	complex_linear grid(complex_linear x);

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
	complex conjugate(complex);

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
	complex sin1(complex);
	complex cos1(complex);

	complex fct(complex);
	complex fctIntegral(complex, complex);
	complex Harmonic(complex);
	complex zeta(complex);

	// infinitesimal (infsim)
	const unsigned int accuracy = 12u;
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
				str += " + ";
				/**/ if (i == acch + 1)		str += "inf";
				else if (i == acch - 1)		str += "0";
				else if (i > acch)	str += "inf^" + std::to_string(i - acch);
				else if (i < acch)	str += "0^" + std::to_string(acch - i);
				else {
					str += "(" + numstr + ")\n";
					continue;
				}
				if (this->Pol[i].R > 0.999999 && this->Pol[i].R < 1.000001 &&
					this->Pol[i].i > -0.000001 && this->Pol[i].i < 0.000001) {
					str += "\n";
					continue;
				}
				str += "*(" + numstr + ")\n";
			}
			if (str == "") return "ConstZero";
			return "   " + str.substr(3);
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
	infsim fct(infsim x);
	infsim getFctIntegralConst();
	infsim fctIntegral(infsim x);
	infsim fctIntegral(infsim x, infsim y);
	infsim Harmonic(infsim x);
	infsim zeta(infsim x);

	typedef complex number;
}