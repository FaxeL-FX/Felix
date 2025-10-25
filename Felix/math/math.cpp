#include "math.h"
namespace math {
	//	number
	const long double
		pi = 3.1415926535897932384626433832795028841971693993751058209749445923,
		inf = 1.0 / 0.0;
	const complex
		i(0, 1),
		fctIntegralConstant = getFctIntegralConstant();
	const infsim
		infinity(acch + 1, 1.0),
		zero(acch - 1, 1.0),
		lnInf = -ln_sum(infinity),
		fctIntegralConst = getFctIntegralConst();

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

	complex getFctIntegralConstant() {
		const int n = 256;
		complex res, a = complex(0, pi / n), b = 8 / 256;
		for (int k = 0; k < n; k++) {
			complex t = exp(k * b - (n - k) * a) + 1;
			res = res + t * exp(-t) / ln(t) * (exp((k + 1) * b - (n - k - 1) * a) - exp(k * b - (n - k) * a));
		}
		return res;
	}

	//	complex
	complex::complex(infsim x) {
		this->R = x.getNum(acch).R;
		this->i = x.getNum(acch).i;
	}

	complex operator+(complex x, complex y) { return complex(x.R + y.R, x.i + y.i); }
	complex operator-(complex x) { return complex(-x.R, -x.i); }
	complex operator-(complex x, complex y) { return complex(x.R - y.R, x.i - y.i); }
	complex operator*(complex x, complex y) { return complex(x.R * y.R - x.i * y.i, x.R * y.i + x.i * y.R); }
	complex operator/(complex x, complex y) {
		if (x == y) return 1;
		if (x == 0) return 0;
		if (0 == y) {
			if (x.i == 0) return complex(x.R / 0.0);
			if (x.R == 0) return complex(0, x.i / 0.0);
			return complex(x.R / 0.0, x.i / 0.0);
		}
		long double denominator = 1.0 / (y.R * y.R + y.i * y.i);
		return complex((x.R * y.R + x.i * y.i) * denominator, (y.R * x.i - x.R * y.i) * denominator);
	}
	complex operator%(complex x, complex y) { return x - y * floor(x / y); }
	complex operator<<(complex x, int n) {
		if (n < 0) {
			for (int i = 0; i < -n / 64; i++) x = x / ((unsigned long long)1 << 63) * 0.5;
			return x / ((unsigned long long)1 << -n % 64);
		}
		for (int i = 0; i < n / 64; i++) x = x * ((unsigned long long)1 << 63) * 2;
		return x * ((unsigned long long)1 << n % 64);
	}																						//	need to change for "flex_float" instead "long double"
	complex operator>>(complex x, int n) { return x << -n; }
	bool operator==(complex x, complex y) {
		if (x.R == y.R && x.i == y.i) return true;
		return false;
	}

	//	functions
	long double rand(int seed, std::vector<long double> v) {
		long double res = cos(seed);
		for (int i = 1; i <= 16; i++)
			for (auto n : v)
				res = res + sin(n) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}
	long double rand(int seed, std::vector<complex> v) {
		long double res = cos(seed);
		for (int i = 1; i <= 16; i++)
			for (auto n : v)
				res = res + sin(((complex)n).R) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}
	long double rand(int seed, complex z) {
		std::vector<long double> v = { z.R, z.i };
		return rand(seed, v);
	}
	long double rand(int n) {
		std::vector<long double> v = { (long double)n };
		return rand(n, v);
	}

	//	long double
	long double floor(long double x) { return std::floor(x); }
	long double ceil(long double x) { return std::ceil(x); }
	long double sign(long double x) {
		if (x < 0) return -1;
		return 1;
	}
	long double round(long double x) {
		return (long long)(x + 0.5);
	}

	long double exp(long double x) {
		long double res = x / ((unsigned long long)1 << 32);
		if (res < 0) res = -res;
		for (int i = 0; i < 32; i++) res = 2 * res + res * res;
		if (x < 0) return 1 / (res + 1);
		return res + 1;
	}
	long double log(long double x) {
		long double res = 0, part = 1;
		for (int k = 1; k < 32; k++) {
			part = part * (1 - x);
			res = res + (1 - part) / (k * ((unsigned long long)1 << k));
		}
		return res;
	}
	long double ln(long double x) { return log(fract(x)) + 0.6931471805599453 * E(x); }
	long double pow(long double x, int n) {
		long double res = 1;
		if (n > 0) for (;0 < n;n--) {
			res *= x;
		}
		if (n < 0) for (; n < 0; n++) {
			res /= x;
		}
		return res;
	}

	long double sqrt(long double x) { return exp(0.5 * ln(x)); }
	long double inv_sqrt(long double x) { return exp(-0.5 * ln(x)); }

	long double cos(long double x) {
		long double res = x / ((unsigned long long)1 << 16);
		res = 2 - res * res;
		for (int i = 0; i < 16; i++) res = res * res - 2;
		return res * 0.5;
	}
	long double sin(long double x) {
		long double res = x, sqr;
		for (int i = 1; i <= 32; i++) {
			sqr = x / (pi * i);
			res *= 1 - sqr * sqr;
		}
		return res * exp(-sqr * sqr * 32);
	}
	long double arccos(long double x) {
		if (x == 1) return 0;
		if (x < 0) return pi - arccos(-x);
		long double res = sqrt((x + 1) * 2);
		for (int i = 1; i < 8; i++) res = sqrt(2 + res);
		return sqrt(2 - res) * (1 << 8);
	}

	double factorial(double x) {
		double res = 1.0;
		for (double k = 2.0; k <= x; k++) res *= k;
		for (double k = -1.0; x < k; k--) res *= k;
		if (x < -1.0) return 1.0 / res;
		return res;
	}
	double Binomial(double n, double k) {
		double res = n;
		for (int l = 1; l <= k; l++) {
			res *= n / l - 1;
		}
		return res;
	}

	//	complex
	complex floor(complex x) { return complex(floor(x.R), floor(x.i)); }
	complex ceil(complex x) { return complex(ceil(x.R), ceil(x.i)); }
	long double abs(complex x) {
		if (x.i == 0) return sign(x.R) * x.R;
		if (x.R == 0) return sign(x.i) * x.i;
		return sqrt(x.R * x.R + x.i * x.i);
	}
	long double inv_abs(complex x) { return inv_sqrt(x.R * x.R + x.i * x.i); }
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
	long double arg(complex x) {
		long double cosine;
		if (x.i == 0)
			if (x.R < 0)	cosine = -1;
			else			cosine = 1;
		else				cosine = normalize(x).R;
		if (x.i < 0) return -arccos(cosine);
		return arccos(cosine);
	}
	bool exist(complex x) {
		if (x.R == x.R && x.i == x.i) return true;
		return false;
	}

	complex exp(complex x) {
		if (x.R < -709) return 0;
		complex res = x >> 64;
		if (res.R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.R < 0) return 1 / (res + 1);
		return res + 1;
	}
	complex ln(complex x) { return complex(ln(abs(x)), arg(x)); }
	complex pow(complex x, complex y) { return exp(y * ln(x)); }
	complex sqrt(complex x) { return exp(0.5 * ln(x)); }
	complex inv_sqrt(complex x) { return exp(-0.5 * ln(x)); }

	complex Re(complex x) { return x.R; }
	complex Im(complex x) { return x.i; }
	complex grid(complex x) {
		complex res = x - math::floor(x + complex(0.5, 0.5));
		if (res.R < 0) res.R = -res.R;
		if (res.i < 0) res.i = -res.i;
		if (res.R > res.i) return (res.i - res.i * res.i) * 4;
		return (res.R - res.R * res.R) * 4;
	}

	complex conjugate(complex x) { return complex(x.R, -x.i); }

	complex cosh(complex x) { return (exp(x) + exp(-x)) * 0.5; }
	complex sinh(complex x) { return (exp(x) - exp(-x)) * 0.5; }
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
	complex arctanh(complex x) { return -0.5 * ln((1 - x) / (1 + x)); }

	complex cos(complex x) { return        cosh(mul_i(x)); }
	complex sin(complex x) { return -mul_i(sinh(mul_i(x))); }
	complex cot(complex x) { return  mul_i(coth(mul_i(x))); }
	complex tan(complex x) { return -mul_i(tanh(mul_i(x))); }
	complex arccos(complex x) { return -mul_i(arccosh(x)); }
	complex arcsin(complex x) { return -mul_i(arcsinh(mul_i(x))); }
	complex arccot(complex x) { return  mul_i(arccoth(mul_i(x))); }
	complex arctan(complex x) { return -mul_i(arctanh(mul_i(x))); }

	complex sin1(complex x) {
		x.R = fmod(x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R > 0.5) return sin1(1 - x);
		else if (x.R < -0.5) return sin1(-1 - x);
		complex res = x;
		for (int i = 1; i < 64; i++) {
			res = res * (1 - x * x / (i * i));
		}
		return res;
	}
	complex cos1(complex x) {
		x.R = fmod(x.R + 1.0, 2.0) - 1.0;
		/**/ if (x.R > 0.5) return -cos1(1 - x);
		else if (x.R < -0.5) return -cos1(-1 - x);
		complex res = 1;
		for (double i = 0.5; i < 64; i++) {
			res = res * (1 - x * x / (i * i));
		}
		return res;
	}
	complex Binom(complex n, complex k) {
		if (k.i == 0 && k.R == floor(k.R)) {
			complex res = 1;
			if (0 <= k.R) {
				for (int l = 0; l < k.R; l++) {
					res = res * (n - l) / (l + 1);
				}
			}
			return res;
		}
		return fct(n) / (fct(k) * fct(n - k));
	}

	complex fct(complex x) {
		if (x.R < -0.5) return -1 / (sin1(x) * fct(-1 - x));
		//if (x.R > 0) return x * fct(x - 1);
		float n = 2048.5; // 2048.5
		complex res = exp(-x + (x + 0.5) * ln(x + n) - 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * i / n * (x + n) / (x + i);
		return res;
	}
	complex fctIntegral(complex x, complex N) {
		int n = 256;
		double invN = 1 / (double)n;
		complex res, t = ln(complex(-1.0)), dt, w = 0.0;
		dt = (ln(n) - t) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			w = exp(t) + 1.0;
			res = res + exp(-w + x * ln(w) - N * ln(ln(w))) * (exp(t + 0.5 * dt) - exp(t - 0.5 * dt));
		}
		return res;

		/*const int n = 256;
		complex res = fctIntegralConstant, logX = ln(x) / n;
		for (int k = 0; k < n; k++) {
			res = res + fct(exp(k * logX)) * (exp((k + 1) * logX) - exp(k * logX));
		}
		return res;*/
	}
	complex Harmonic(complex x) {
		if (x.R < -0.5) return Harmonic(-x - 1) - cos1(x) / sin1(x);
		int n = 2048;
		complex res = ln(1 + x / n);
		for (int k = 1; k <= n; k++)
			res = res + x / (k * (x + k));
		return res;
	}
	complex zeta(complex x) {
		if (x.R > 3) return -fct(x) * exp(-x * 1.83787706640934548356) * sin1(x * 0.5) * zeta(-x - 1);
		if (x.R < -170) return 1;
		int n = 64;
		long double pow2 = 1;
		complex res1 = 0;

		for (int i = 0; i < n; i++) {
			complex res2 = 0;
			for (int j = 0; j <= i; j++) {
				res2 = res2 + (1 - 2 * (j % 2)) * factorial(i) / factorial(i - j) * pow(j + 1, x) / factorial(j);
			}
			pow2 *= 0.5;
			res1 = res1 + res2 * pow2;
		}

		return res1 / (1 - 2 * pow(2, x));
	}
	complex zetaByFct(complex x) {
		if (x == complex(-1)) return -1;
		if (x.R < -1) {
			complex res = 0;
			for (int k = 1; k < 64; k++) {
				res = res + (2 * (k % 2) - 1) * pow(k, x);
			}
			return -res / (1 - 2 * pow(2, x)) * sin1(x) * fct(-x - 1);
		}

		complex res = 0;
		for (int k = 1; k < x.R + 64; k++) {
			res = res + zetaByFct(x - k) * (2 * (k % 2) - 1) / factorial(k + 1);
		}
		return res;
	}

	// infsim
	infsim operator+(infsim x, infsim y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) + y.getNum(i));
		return x;
	}
	infsim operator-(infsim x) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, -x.getNum(i));
		return x;
	}
	infsim operator-(infsim x, infsim y) { return x + -y; }
	infsim operator*(infsim x, infsim y) {
		infsim res;
		for (int i = 0; i < accuracy; i++)
			for (int j = acch - i; i + j - acch < accuracy; j++)
				res.setNum(i + j - acch, res.getNum(i + j - acch) + x.getNum(i) * y.getNum(j));
		return res;
	}
	infsim inv(infsim x) {
		int maxPower = accuracy - 1;
		for (; 0 <= maxPower && x.getNum(maxPower).isZero(); maxPower--) {}
		int invPower = acch - maxPower;
		complex denominator = 1.0 / (x.getNum(maxPower));
		x =  mul(x * infsim(acch + invPower, 1.0), denominator);
		infsim numerator(acch + invPower, 1.0), res;
		for (int i = acch + invPower; 0 <= i; i--) {
			infsim mul(i, numerator.getNum(i));
			res = res + mul;
			numerator = numerator - mul * x;
		}
		return mul(res, denominator);
	}
	infsim operator/(infsim x, infsim y) { return x * inv(y); }
	infsim mul(infsim x, complex y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) * y);
		return x;
	}
	infsim div(infsim x, complex y) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) / y);
		return x;
	}
	infsim operator<<(infsim x, int n) {
		if (n < 0)	for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) >> -n);
		else		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) << n);
		return x;
	}
	infsim operator>>(infsim x, int n) { return x << -n; }


	infsim Re(infsim x) {
		infsim res;
		for (int i = 0; i < accuracy; i++) res.setNum(i, x.getNum(i).R);
		return res;
	}
	infsim Im(infsim x) {
		infsim res;
		for (int i = 0; i < accuracy; i++) res.setNum(i, x.getNum(i).i);
		return res;
	}
	complex grid(infsim x) {
		complex res = x.getNum(acch) - math::floor(x.getNum(acch) + complex(0.5, 0.5));
		if (res.R < 0) res.R = -res.R;
		if (res.i < 0) res.i = -res.i;
		if (res.R > res.i) return res.i;
		return res.R;
	}
	infsim mul_i(infsim x) {
		for (int i = 0; i < accuracy; i++) x.setNum(i, x.getNum(i) * math::i);
		return x;
	}
	infsim floor(infsim x) {
		x.Pol[acch] = floor(x.Pol[acch]);
		return x;
	}
	infsim ceil(infsim x) {
		x.Pol[acch] = ceil(x.Pol[acch]);
		return x;
	}

	infsim exp(infsim x) {
		infsim res = x >> 64;
		if (res.getNum(acch).R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.getNum(acch).R < 0) return 1 / (res + 1);
		return res + 1;
	}
	infsim ln_sum(infsim x) {
		complex c = x.getNum(acch);
		infsim res, s = 1, t = 1 - infsim(c) / x;
		for (int k = 1; k < 32; k++) {
			s = s * t;
			res = res + div(s, k);
		}
		return res + infsim(ln(c));
	}
	infsim ln_integral(infsim x) {
		int n = 128;
		double invN = 1 / (double)n;
		infsim res, t = 1.0, th = i * abs(x.getNum(acch)), dt;
		if (x.getNum(acch).i < 0.0) th = -th;
		if (x.getNum(acch).R < 0.0) {
			dt = (th - 1.0) * invN;
			for (int i = 0; i < n; i++) {
				t = t + dt;
				res = res + inv(t) * dt;
			}
		}
		else th = 1.0;
		dt = (x - th) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			res = res + inv(t) * dt;
		}
		return res;
	}
	infsim ln(infsim x) {
		int maxPower = accuracy - 1;
		for (maxPower; 0 <= maxPower && x.getNum(maxPower).isZero(); maxPower--) {}
		int maxPowerVal = maxPower - acch;
		x = infsim(acch - maxPowerVal, 1.0) * x;
		//return mul(lnInf, maxPowerVal) + ln_integral(x);
		return mul(lnInf, maxPowerVal) + ln_sum(x);
	}
	infsim pow(infsim x, infsim y) {
		if (y.getNum(acch).R >= 1) return pow(x, y - 1) * x;
		if (y.getNum(acch).R < 0) return pow(x, y + 1) / x;
		return exp(y * ln(x));
	}
	infsim sqrt(infsim x) { return exp(0.5 * ln(x)); }
	infsim inv_sqrt(infsim x) { return exp(-0.5 * ln(x)); }

	infsim cosh(infsim x) { return (exp(x) + exp(-x)) * 0.5; }
	infsim sinh(infsim x) { return (exp(x) - exp(-x)) * 0.5; }
	infsim coth(infsim x) {
		infsim positive = exp(x), negative = exp(-x);
		return (positive + negative) / (positive - negative);
	}
	infsim tanh(infsim x) {
		infsim positive = exp(x), negative = exp(-x);
		return (positive - negative) / (positive + negative);
	}
	infsim arccosh(infsim x) { return ln(x + sqrt(x * x - 1)); }
	infsim arcsinh(infsim x) { return ln(x + sqrt(x * x + 1)); }
	infsim arccoth(infsim x) { return arctanh(1 / x); }
	infsim arctanh(infsim x) { return -0.5 * ln((1 - x) / (1 + x)); }

	infsim cos(infsim x) { return      cosh(mul_i(x)); }
	infsim sin(infsim x) { return -mul_i(sinh(mul_i(x))); }
	infsim cot(infsim x) { return  mul_i(coth(mul_i(x))); }
	infsim tan(infsim x) { return -mul_i(tanh(mul_i(x))); }
	infsim arccos(infsim x) { return -mul_i(arccosh(x)); }
	infsim arcsin(infsim x) { return -mul_i(arcsinh(mul_i(x))); }
	infsim arccot(infsim x) { return  mul_i(arccoth(mul_i(x))); }
	infsim arctan(infsim x) { return -mul_i(arctanh(mul_i(x))); }

	infsim sin1(infsim x) {
		x.setNum(acch, fmod(x.getNum(acch).R + 1.0, 2.0) - 1.0);
		/**/ if (x.getNum(acch).R > 0.5) return sin1(1 - x);
		else if (x.getNum(acch).R < -0.5) return sin1(-1 - x);
		infsim res = x;
		for (int i = 1; i < 64; i++) {
			res = res * (1 - div(x * x, i * i));
		}
		return res;
	}
	infsim cos1(infsim x) {
		x.setNum(acch, fmod(x.getNum(acch).R + 1.0, 2.0) - 1.0);
		/**/ if (x.getNum(acch).R > 0.5) return -cos1(1 - x);
		else if (x.getNum(acch).R < -0.5) return -cos1(-1 - x);
		infsim res = 1;
		for (double i = 0.5; i < 64; i++) {
			res = res * (1 - div(x * x, i * i));
		}
		return res;
	}
	infsim Binom(infsim n, infsim k) {
		return fct(n) / (fct(k) * fct(n - k));
	}

	infsim fct(infsim x) {
		if (x.getNum(acch).R <= -1) return fct(x + 1) / (x + 1);
		if (x.getNum(acch).R > 0) return fct(x - 1) * x;
		float n = 512;
		infsim res = exp(-x + (x + 0.5) * ln(x + n) - 0.5 * ln(n));
		for (int i = 1; i < n; i++) res = res * i / n * (x + n) / (x + i);
		return res;
	}
	infsim getFctIntegralConst() {
		const int n = 256;
		infsim res, a = infsim(0, pi / n), b = 8 / 256;
		for (int k = 0; k < n; k++) {
			infsim t = exp(k * b - (n - k) * a) + 1;
			res = res + t * exp(-t) / ln(t) * (exp((k + 1) * b - (n - k - 1) * a) - exp(k * b - (n - k) * a));
		}
		return res;
	}
	infsim fctIntegral(infsim x) {
		const int n = 16;
		infsim res = fctIntegralConstant, logX = ln(x) / n;
		for (int k = 0; k < n; k++) {
			res = res + fct(exp(k * logX)) * (exp((k + 1) * logX) - exp(k * logX));
		}
		return res;
	}
	infsim fctIntegral(infsim x, infsim y) {
		int n = 32;
		double invN = 1 / (double)n;
		infsim res, t = ln(infsim(-1.0)), dt, w = 0.0;
		dt = (infsim(256) - t) * invN;
		for (int i = 0; i < n; i++) {
			t = t + dt;
			w = exp(t) + 1.0;
			res = res + exp(-w + x * ln(w) - y * ln(ln(w))) * (exp(t) - exp(t - dt));
		}
		return res;
	}
	infsim Harmonic(infsim x) {
		if (x.getNum(acch).R < -0.5) return Harmonic(-x - 1) - pi * cot(pi * x);
		int n = 64;
		infsim res = ln(1 + x / n);
		for (int k = 1; k <= n; k++)
			res = res + x / (k * (x + k));
		return res;
	}
	infsim zeta(infsim x) {
		if (x.getNum(acch).R > 0) return -pow(2 * pi, -x) * sin1(x * 0.5) * fct(x) * zeta(-x - 1);
		infsim res1 = 0;
		for (int i = 0; i < 12; i++) {
			infsim res2 = 0;
			for (int j = 0; j <= i; j++) {
				res2 = res2 + (1 - 2 * (j % 2)) * factorial(i) * pow(j + 1, x) / (factorial(j) * factorial(i - j));
			}
			res1 = res1 + mul(res2, pow(complex(2.0), -(i + 1)));
		}
		return res1 / (1 - pow(2, x + 1));
	}
}