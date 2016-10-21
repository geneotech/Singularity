#pragma once
extern int is_prime_lookup[104730];
void is_prime_lookup_init();
bool is_prime_lookup_safe(int num);
int nth_prime(int n);
int max_primes();