#pragma once
#include <vector>
#include <algorithm>
#include <cassert>
#include "primes.h"

typedef std::vector<int> factor_vector;
bool operator<(const factor_vector& a, const factor_vector& b);
typedef unsigned long long ull;

struct full_factors_info {
	ull number;

	int total_factors = 0;
	std::vector<int> prime_factors_sequence;
	std::vector<int> prime_ordinals_sequence;
	std::vector<int> prime_factors;
	std::vector<int> prime_ordinals;
	std::vector<int> exponents;

	std::vector<int> list_consecutive_exponents();
	void push_factor(int new_prime_factor);
	bool prime();

	bool semiprime();

	bool twin_prime();

	bool prime_power();

	bool prime_power_gt_2();

	bool prime_square();

	bool has_unique_twofactorization();

	int smaller_twofactor();

	int bigger_twofactor();
};

std::pair<int, int> closest_divisor_pair(int n);
full_factors_info get_full_factor_info(ull n);