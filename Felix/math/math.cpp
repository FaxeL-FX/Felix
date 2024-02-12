#include "math.h"
namespace math {
	//	number
	const long double
		pi = 3.1415926535897932384626433832795028841971693993751058209749445923,
		inf = 1.0 / 0.0;

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

	fixed_point32 operator+(fixed_point32 x, fixed_point32 y) {
		x.num = x.num + y.num;
		return x;
	}
	fixed_point32 operator-(fixed_point32 x) {
		x.num = -x.num;
		return x;
	}
	fixed_point32 operator-(fixed_point32 x, fixed_point32 y) { return x + -y; }
	fixed_point32 operator*(fixed_point32 x, fixed_point32 y) {
		x.num = ((x.num >> 32) * y.num) << 32 + (((x.num << 32) >> 32) * y.num) >> 32;
		return x;
	}
	fixed_point32 operator<<(fixed_point32 x, int n) {
		if (n < 0)  x.num = x.num >> n;
		else		x.num = x.num << n;
		return x;
	}
	fixed_point32 operator>>(fixed_point32 x, int n) { return x << -n; }
	bool operator==(fixed_point32 x, fixed_point32 y) {
		if (x.num == y.num) return true;
		return false;
	}

	complex_linear getFctIntegralConstant() {
		const int n = 256;
		complex_linear res, a = complex_linear(0, pi / n), b = 8 / 256;
		for (int k = 0; k < n; k++) {
			complex_linear t = exp(k * b - (n - k) * a) + 1;
			res = res + t * exp(-t) / ln(t) * (exp((k + 1) * b - (n - k - 1) * a) - exp(k * b - (n - k) * a));
		}
		return res;
	}

	//	complex
	complex_linear::complex_linear(complex_exponential x) {
		this->R = x.r * cos(x.a);
		this->i = x.r * sin(x.a);
	}
	complex_exponential::complex_exponential(long double x) {
		this->r = x * sign(x);
		if (x >= 0)	this->a = 0;
		else		this->a = pi;
	}
	complex_exponential::complex_exponential(complex_linear x) {
		this->r = abs(x);
		this->a = arg(x);
	}
	const complex_linear
		i(0, 1),
		fctIntegralConstant = getFctIntegralConstant();

	complex_linear operator+(complex_linear x, complex_linear y) { return complex_linear(x.R + y.R, x.i + y.i); }
	complex_linear operator-(complex_linear x) { return complex_linear(-x.R, -x.i); }
	complex_linear operator-(complex_linear x, complex_linear y) { return complex_linear(x.R - y.R, x.i - y.i); }
	complex_linear operator*(complex_linear x, complex_linear y) { return complex_linear(x.R * y.R - x.i * y.i, x.R * y.i + x.i * y.R); }
	complex_linear operator/(complex_linear x, complex_linear y) {
		if (x == y) return 1;
		if (x == 0) return 0;
		if (0 == y) {
			if (x.i == 0) return complex_linear(x.R / 0.0);
			if (x.R == 0) return complex_linear(0, x.i / 0.0);
			return complex_linear(x.R / 0.0, x.i / 0.0);
		}
		long double denominator = 1.0 / (y.R * y.R + y.i * y.i);
		return complex_linear((x.R * y.R + x.i * y.i) * denominator, (y.R * x.i - x.R * y.i) * denominator);
	}
	complex_linear operator%(complex_linear x, complex_linear y) { return x - y * floor(x / y); }
	complex_linear operator<<(complex_linear x, int n) {
		if (n < 0) {
			for (int i = 0; i < -n / 64; i++) x = x / ((unsigned long long)1 << 63) * 0.5;
			return x / ((unsigned long long)1 << -n % 64);
		}
		for (int i = 0; i < n / 64; i++) x = x * ((unsigned long long)1 << 63) * 2;
		return x * ((unsigned long long)1 << n % 64);
	}																						//	need to change for "flex_float" instead "long double"
	complex_linear operator>>(complex_linear x, int n) { return x << -n; }
	bool operator==(complex_linear x, complex_linear y) {
		if (x.R == y.R && x.i == y.i) return true;
		return false;
	}

	complex_exponential operator+(complex_exponential x, complex_exponential y) {
		long double
			radius = sqrt(x.r * x.r + y.r * y.r + 2 * x.r * y.r * cos(y.a - x.a)),
			angle = sign(x.r * sin(x.a) + y.r * sin(y.a)) * arccos((x.r * cos(x.a) + y.r * cos(y.a)) / radius);
		return complex_exponential(radius, angle);
	}
	complex_exponential operator-(complex_exponential x) { return complex_exponential(x.r, x.a + pi); }
	complex_exponential operator-(complex_exponential x, complex_exponential y) { return x + -y; }
	complex_exponential operator*(complex_exponential x, complex_exponential y) { return complex_exponential(x.r * y.r, x.a + y.a); }
	complex_exponential operator/(complex_exponential x, complex_exponential y) { return complex_exponential(x.r / y.r, x.a - y.a); }
	complex_exponential operator%(complex_exponential x, complex_exponential y) { return x - y * floor(x / y); }
	complex_exponential operator<<(complex_exponential x, int n) {
		if (n < 0) {
			for (int i = 0; i < -n / 64; i++) x = x / ((unsigned long long)1 << 63) * 0.5;
			return x / ((unsigned long long)1 << -n % 64);
		}
		for (int i = 0; i < n / 64; i++) x = x * ((unsigned long long)1 << 63) * 2;
		return x * ((unsigned long long)1 << n % 64);
	}
	complex_exponential operator>>(complex_exponential x, int n) { return x << -n; }
	bool operator==(complex_exponential x, complex_exponential y) {
		if (x.r == y.r && x.a == y.a) return true;
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
				res = res + sin(((complex_linear)n).R) + res - 16 * floor(res * 0.0625);
		res = fmod(res, 1.0);
		return res;
	}
	long double rand(int seed, complex_linear z) {
		std::vector<long double> v = { z.R, z.i };
		return rand(seed, v);
	}
	long double rand(int n) {
		std::vector<long double> v = { (long double)n };
		return rand(n, v);
	}

	//	long double
	long double floor(long double x) {
		return std::floor(x);
	}
	long double ceil(long double x) {
		return std::ceil(x);
	}
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

	//	complex_linear
	std::string toString(complex_linear x) {
		if (x.i == 0) x.i = 0;
		if (x.i < 0) return std::to_string(x.R) + std::to_string(x.i) + "i";
		return std::to_string(x.R) + "+" + std::to_string(x.i) + "i";
	}

	complex_linear floor(complex_linear x) { return complex_linear(floor(x.R), floor(x.i)); }
	complex_linear ceil(complex_linear x) { return complex_linear(ceil(x.R), ceil(x.i)); }
	long double abs(complex_linear x) {
		if (x.i == 0) return sign(x.R) * x.R;
		if (x.R == 0) return sign(x.i) * x.i;
		return sqrt(x.R * x.R + x.i * x.i);
	}
	long double inv_abs(complex_linear x) { return inv_sqrt(x.R * x.R + x.i * x.i); }
	complex_linear normalize(complex_linear x) {
		if (x.i == 0)
			if (x.R < 0)	return -1;
			else			return 1;
		if (x.R == 0)
			if (x.i < 0)	return -i;
			else			return i;
		return x * inv_sqrt(x.R * x.R + x.i * x.i);
	}
	complex_linear mul_i(complex_linear x) { return complex_linear(-x.i, x.R); }
	long double arg(complex_linear x) {
		long double cosine;
		if (x.i == 0)
			if (x.R < 0)	cosine = -1;
			else			cosine = 1;
		else				cosine = normalize(x).R;
		if (x.i < 0) return -arccos(cosine);
		return arccos(cosine);
	}
	bool exist(complex_linear x) {
		if (x.R == x.R && x.i == x.i) return true;
		return false;
	}

	complex_linear exp(complex_linear x) {
		if (x.R < -709) return 0;
		complex_linear res = x >> 64;
		if (res.R < 0) res = -res;
		for (int i = 0; i < 64; i++) res = (res << 1) + res * res;
		if (x.R < 0) return 1 / (res + 1);
		return res + 1;
	}
	complex_linear ln(complex_linear x) { return complex_linear(ln(abs(x)), arg(x)); }
	complex_linear pow(complex_linear x, complex_linear y) {
		if (y == 0) return 1;
		if (y.R < 0)	return 1 / pow(x, -y);
		return exp(y * ln(x));
	}
	complex_linear sqrt(complex_linear x) { return exp(0.5 * ln(x)); }
	complex_linear inv_sqrt(complex_linear x) { return exp(-0.5 * ln(x)); }

	complex_linear fct(complex_linear x) {
		if (x.R < -0.5) return -pi / ((complex_linear)sin(pi * x) * fct(-1 - x));
		float n = 32.5;
		complex_linear res = exp(-x + (x + n) * ln(x + n) - n * ln(n));
		for (int i = 1; i < n; i++) res = res * (i / (x + i));
		return res;
	}

	//	complex_exponential
	std::string toString(complex_exponential x) { return "[" + std::to_string(x.r) + ":" + std::to_string(x.a) + "]"; }
	complex_exponential floor(complex_exponential x) { return floor((complex_linear)x); }
	long double abs(complex_exponential x) { return x.r; }
	long double inv_abs(complex_exponential x) { return 1 / x.r; }
	complex_exponential normalize(complex_exponential x) { return complex_exponential(1, x.a); }
	complex_exponential mul_i(complex_exponential x) { return complex_exponential(x.r, x.a + (0.5 * pi)); }
	long double arg(complex_exponential x) { return x.a; }
	bool exist(complex_exponential x) {
		if (x.r == x.r) return true;
		return false;
	}

	complex_exponential exp(complex_exponential x) { return complex_exponential(exp(x.r * cos(x.a)), x.r * sin(x.a)); }
	complex_exponential ln(complex_exponential x) { return complex_linear(ln(x.r), x.a); }
	complex_exponential pow(complex_exponential x, complex_exponential y) {
		long double c = cos(y.a), s = sin(y.a), l = ln(x.r);
		return complex_exponential(exp(y.r * (l * c - (long double)x.a * s)), y.r * (l * s + (long double)x.a * c));
	}
	complex_exponential sqrt(complex_exponential x) { return complex_exponential(sqrt(x.r), x.a * 0.5); }
	complex_exponential inv_sqrt(complex_exponential x) { return complex_exponential(inv_sqrt(x.r), x.a * (-0.5)); }

	complex_exponential fct(complex_exponential x) {
		complex_linear X = x, productRes = 1;
		if (X.R < -0.5) return -pi / ((complex_exponential)sin(pi * x) * fct((complex_exponential)(-1 - X)));
		float n = 32.5;
		complex_exponential Xn = X + n, res = pow(Xn, Xn) / exp(x) * exp(-n * ln(n));
		for (int i = 1; i < n; i++) productRes = productRes * i  / (X + i);
		return res * (complex_exponential)productRes;
	}

	//	complex
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

	complex fctIntegral(complex x, complex N) {
		const int n = 256;
		complex res = fctIntegralConstant, logX = ln(x) / n;
		for (int k = 0; k < n; k++) {
			res = res + fct(exp(k * logX)) * (exp((k + 1) * logX) - exp(k * logX));
		}
		return res;
	}
	complex Harmonic(complex x) {
		if (x.R < -0.5) return Harmonic(-x - 1) - pi * cot(pi * x);
		int n = 64;
		complex res = ln(1 + x / n);
		for (int k = 1; k <= n; k++)
			res = res + x / (k * (x + k));
		return res;
	}
}