#pragma once
#include <vector>
#include <algorithm>
#include <cassert>
#include "primes.h"

typedef std::vector<int> factor_vector;
bool operator<(const factor_vector& a, const factor_vector& b) {
	return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

typedef unsigned long long ull;

struct full_factors_info {
	ull number;

	int total_factors = 0;
	std::vector<int> prime_factors_sequence;
	std::vector<int> prime_ordinals_sequence;
	std::vector<int> prime_factors;
	std::vector<int> prime_ordinals;
	std::vector<int> exponents;

	std::vector<int> list_consecutive_exponents() {
		int length = *prime_ordinals.rbegin();

		std::vector<int> out;
		out.resize(length);

		std::fill(out.begin(), out.end(), 0);

		for (int i = 0; i < prime_ordinals.size(); ++i) {
			out[prime_ordinals[i]-1] = exponents[i];
		}

		return out;
	}

	void push_factor(int new_prime_factor) {
		if (prime_factors.empty() || *prime_factors.rbegin() != new_prime_factor) {
			prime_ordinals.push_back(new_prime_factor < max_primes() ? is_prime_lookup[new_prime_factor] : 0);
			prime_factors.push_back(new_prime_factor);
			exponents.push_back(1);
		}
		else
		{
			(*exponents.rbegin())++;
		}

		prime_factors_sequence.push_back(new_prime_factor);
		prime_ordinals_sequence.push_back(new_prime_factor < max_primes() ? is_prime_lookup[new_prime_factor] : 0);
		++total_factors;
	}

	bool prime() {
		return prime_factors_sequence.size() == 1;
	}

	bool semiprime() {
		return prime_factors.size() == 2 && prime_factors_sequence.size() == 2;
	}

	bool prime_power() {
		return prime_factors.size() == 1;
	}

	bool prime_power_gt_2() {
		return prime_factors.size() == 1 && exponents[0] > 1;
	}

	bool prime_square() {
		return prime_factors.size() == 1 && exponents[0] == 2;
	}

	bool has_unique_twofactorization() {
		auto& f = prime_factors_sequence;

		return f.size() == 2 || (f.size() == 3 && f[0] == f[1] && f[1] == f[2]);
	}

	int smaller_twofactor() {
		auto& f = prime_factors_sequence;

		if (f.size() == 2) {
			return f[0];
		}
		else if (f.size() == 3) {
			return f[0];
		}
		else assert(0);
	}

	int bigger_twofactor() {
		auto& f = prime_factors_sequence;

		if (f.size() == 2) {
			return f[1];
		}
		else if (f.size() == 3) {
			return f[0]*f[1];
		}
		else assert(0);
	}
};

std::pair<int, int> closest_divisor_pair(int n) {
	int testNum = (int)sqrt(n);
	while (n % testNum != 0) {
		testNum--;
	}

	return{ testNum, n / testNum };
}

full_factors_info get_full_factor_info(ull n) {
	int z = 2;

	full_factors_info out;
	out.number = n;

	if (n == 1)
	{
		out.push_factor(1);
		return out;
	}

	int upper_limit = (int)sqrt(n) + 1;

	while (z <= upper_limit)
	{
		if (n % z == 0)
		{
			out.push_factor(z);
			n /= z;
		}
		else
			++z;
	}

	if (n > 1)
	{
		out.push_factor(n);
	}

	return out;
}
