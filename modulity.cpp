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
gcd_length(long long a, long long b)
{
	long long l = 0;
	long long c;
	while (a != 0) {
		c = a; a = b%a;  b = c;
		l++;
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