#include "modulity.h"

int modulity(int n, int d) {
	int moduli_cnt = 0;

	while (d) {
		auto res = n % d;

		if (!res)
			break;

		d = res;
		++moduli_cnt;
	};

	return moduli_cnt;
}

int
gcd_length(int a, int b)
{
	int l = 0;
	int c;
	while (a != 0) {
		c = a; a = b%a;  b = c;
		l++;
	}
	return l;
}