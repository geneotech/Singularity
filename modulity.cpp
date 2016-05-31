#include "modulity.h"

long long modulity(long long n, long long d) {
	long long moduli_cnt = 0;

	while (d) {
		auto res = n % d;

		if (!res)
			break;

		d = res;
		++moduli_cnt;
	};

	return moduli_cnt;
}

long long
gcd_length(long long u, long long v)
{
	unsigned long long l = 0;
	while (v != 0) {
		long long r = u % v;
		u = v;
		v = r;
		++l;
	}
	return l;
}

long long gcd(long long u, long long v) {
	while (v != 0) {
		long long r = u % v;
		u = v;
		v = r;
	}
	return u;
}

#include <cmath>
double gcdd(double u, double v) {
	while (v != 0) {
		double r = fmod(u, v);
		u = v;
		v = r;
	}
	return u;
}