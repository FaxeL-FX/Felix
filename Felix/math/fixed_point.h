#pragma once
#include "../include.h"
#include <cstddef>
#include <cstdint>

namespace math {
	template<size_t bit>
	struct fixed_point {
		uint64_t mantissa = 0;
		bool sign = 0;

		fixed_point() {}
		fixed_point(long double x) {
			this->sign = x < 0.0;
			if (this->sign) x = -x;
			for (int i = 0; i < bit; i++) x *= 2.0;
			x = std::fmod(x, 18446744073709551616.0);
			this->mantissa = x;
		}
		fixed_point(double x) {
			this->sign = x < 0.0;
			if (this->sign) x = -x;
			for (int i = 0; i < bit; i++) x *= 2.0;
			x = std::fmod(x, 18446744073709551616.0);
			this->mantissa = x;
		}
		fixed_point(long long x) {
			this->sign = x < 0;
			if (this->sign) x = -x;
			this->mantissa = x << bit;
		}
		fixed_point(int x) {
			this->sign = x < 0;
			if (this->sign) x = -x;
			this->mantissa = (long long)x << bit;
		}

		operator long double(){
			return toDouble();
		}

		fixed_point(uint64_t m, bool s) {
			this->mantissa = m;
			this->sign = s;
		}

		template<size_t new_bit>
		fixed_point<new_bit> convert() {
			fixed_point<new_bit> a;
			a.sign = this->sign;
			if (new_bit < bit)	a.mantissa = this->mantissa >> (bit - new_bit);
			else				a.mantissa = this->mantissa << (new_bit - bit);
			return a;
		}

		long double toDouble() {
			long double a = this->mantissa;
			for (int i = 0; i < bit; i++) a /= 2.0;
			if (this->sign) a *= -1.0;
			return a;
		}
		std::string toString(size_t digits) {
			return std::format("{:.{}f}", this->toDouble(), digits);
		}
	};

	template<size_t bit>
	fixed_point<bit> operator+(fixed_point<bit> a, fixed_point<bit> b) {
		if (a.sign == b.sign)				a.mantissa = a.mantissa + b.mantissa;
		else if (a.mantissa < b.mantissa) {
			a.mantissa = b.mantissa - a.mantissa;
			a.sign = b.sign;
		}
		else								a.mantissa = a.mantissa - b.mantissa;
		return a;
	}
	template<size_t bit>
	fixed_point<bit> operator-(fixed_point<bit> a) {
		a.sign = !a.sign;
		return a;
	}
	template<size_t bit>
	fixed_point<bit> operator-(fixed_point<bit> a, fixed_point<bit> b) { return a + -b; }

	template<size_t bit>
	fixed_point<bit> operator<<(fixed_point<bit> a, size_t n) {
		a.mantissa = a.mantissa << n;
		return a;
	}
	template<size_t bit>
	fixed_point<bit> operator>>(fixed_point<bit> a, size_t n) {
		a.mantissa = a.mantissa >> n;
		return a;
	}

	template<size_t bit>
	fixed_point<bit> operator*(fixed_point<bit> a, fixed_point<bit> b) {
		fixed_point<bit> res;
		if (a.sign == b.sign)	res.sign = 0;
		else					res.sign = 1;
		for (int i = 0; i < 64; i++) {
			if ((a.mantissa << (63 - i)) >> 63) {
				int index = i - bit;
				if (index < 0)  res.mantissa = res.mantissa + (b.mantissa >> (-index));
				else			res.mantissa = res.mantissa + (b.mantissa << (index));
			}
		}
		return res;
	}
	template<size_t bit>
	fixed_point<bit> operator/(fixed_point<bit> a, fixed_point<bit> b) {
		fixed_point<bit> res;
		if (a.sign == b.sign)	res.sign = 0;
		else					res.sign = 1;
		int i_max = 0;
		for (; i_max < 64; i_max++) {
			auto c = b.mantissa << 1;
			if (a.mantissa >= c) b.mantissa = c;
			else break;
		}
		for (int i = i_max + bit; 0 <= i; i--) {
			if (a.mantissa >= b.mantissa) {
				a.mantissa = a.mantissa - b.mantissa;
				res.mantissa = res.mantissa + ((uint64_t)1 << i);
			}
			b.mantissa = b.mantissa >> 1;
		}
		return res;
	}

	template<size_t bit>
	bool operator==(fixed_point<bit> a, fixed_point<bit> b) { return a.sign == b.sign && a.mantissa == b.mantissa; }
	template<size_t bit>
	bool operator>(fixed_point<bit> a, fixed_point<bit> b) {
		if (a.sign == b.sign) {
			if (a.sign) return a.mantissa < b.mantissa;
			else		return a.mantissa > b.mantissa;
		}
		else return b.sign;
	}
	template<size_t bit>
	bool operator<(fixed_point<bit> a, fixed_point<bit> b) { return b > a; }
	template<size_t bit>
	bool operator>=(fixed_point<bit> a, fixed_point<bit> b) { return !(a < b); }
	template<size_t bit>
	bool operator<=(fixed_point<bit> a, fixed_point<bit> b) { return !(a > b); }

	template<size_t bit>
	fixed_point<bit> exp(fixed_point<bit> x) {
		if (x.sign) return 1.0 / exp(-x);
		fixed_point<bit> res = 1.0, part = 1.0;
		for (int k = 1; k < 32; k++) {
			part = part * (x / k);
			res = res + part;
		}
		return res;
	}
}