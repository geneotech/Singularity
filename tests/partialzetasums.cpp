// experimental partial zeta sums with consideration for primes and divisors overall
#include "tests.h"
#include "../img.h"
#include "../modulity.h"

#include "../util.h"
#include "../factorizator.h"
const int FMAX = 1000;

// Create an array for memoization
int f[FMAX] = { 0 };

// Returns n'th fuibonacci number using table f[]
int fib(int n)
{
	// Base cases
	if (n == 0)
		return 0;
	if (n == 1 || n == 2)
		return (f[n] = 1);

	// If fib(n) is already computed
	if (f[n])
		return f[n];

	int k = (n & 1) ? (n + 1) / 2 : n / 2;

	// Applyting above formula [Note value n&1 is 1
	// if n is odd, else 0.
	f[n] = (n & 1) ? (fib(k)*fib(k) + fib(k - 1)*fib(k - 1))
		: (2 * fib(k - 1) + fib(k))*fib(k);

	return f[n];
}

inline float nth_triangular(float n) {
	return  n *(n + 1) / 2;
}

#include <numeric>

unsigned long long iter_factorial(unsigned long long n)
{
	unsigned long long ret = 1;
	for (unsigned long long i = 1; i <= n; ++i)
		ret *= i;
	return ret;
}

unsigned long long iter_nn(unsigned long long n)
{
	unsigned long long ret = n;
	for (unsigned long long i = 1; i <= n; ++i)
		ret *= n;
	return ret;
}

inline double nth_length(float n) {
	return
		1.f / pow(2, n);
	//1.f / (sqrt(n));
}

#define xcoord(xxx) (xoff+(xxx)*scale + width / 2)

void partialzetasumsintegers() {

	//int xoff = -width / 1.5;
	//float scale = 12000;

	int offset = 1;
	int width = 1980;
	int height = 1080;
	int xoff = -width/3;

	int samples = 20;
	int imgi = 0;
	float maxsample = 28.f;
	float scale = 1500;
	int granularity = 1;

	std::vector<vec2> traj;
	std::vector<int> traj_meta;
	std::vector<vec2> primal;

	bool outputiter = false;

	for (int yyy = 0; yyy < 13600000; ++yyy) {
		//int y = nth_prime(yyy);
		//
		int y = yyy;

		double yy = y / double(granularity);

		if (y % 20000 == 0) {
			outputiter = true;

			setwh(width, height);

			gprint(20, 20, "i=%x", yy);
		}
		else {
			outputiter = false;
		}

		vec2 prev;
		double total_rotation = 0;

		for (int s = 1; s <= samples; ++s) {
			double rotational_contribution = yy / s * 360;
			total_rotation += rotational_contribution;

			vec2 offset = vec2(1, 0).rotate(rotational_contribution, vec2()).set_length(nth_length(s));

			//if (!(s % 2)) {
			//	offset = -offset;
			//}

			vec2 newpos = prev + offset;

			//if(outputiter)
			//Line(xcoord(prev.x), prev.y*scale + height / 2, xcoord(newpos.x), newpos.y*scale + height / 2, 255);

			prev = newpos;
		}

		auto ffi = get_full_factor_info(y / granularity);
		int num_factors = ffi.prime_factors_sequence.size();

		if (y % granularity == 0 && num_factors == 1) {
			primal.push_back(prev);
		}

		traj.push_back(prev);
		traj_meta.push_back(num_factors);

		prev = traj[0];

		if (outputiter) {
			for (int i = 0; i < traj.size(); ++i) {
				auto t = traj[i];
				//Line(prev.x*scale + width / 2, prev.y*scale + height / 2, t.x*scale + width / 2, t.y*scale + height / 2, 0, 255, 0);
				setpix(xcoord(t.x), t.y*scale + height / 2, 255.f/traj_meta[i]);
				prev = t;
			}

			if (primal.size() > 1) {
				prev = primal[0];

				for (size_t i = 1; i < primal.size(); ++i) {
					auto p = primal[i];
					//Line(prev.x*scale + width / 2, prev.y*scale + height / 2, p.x*scale + width / 2, p.y*scale + height / 2, 255, 0, 0);
					setpix(xcoord(p.x), p.y*scale + height / 2, 255, 0, 0);
					prev = p;
				}
			}

			makecross(xcoord(0), height / 2, 0, 0, 255);
			makecross(xcoord(1), height / 2, 0, 0, 255);

			std::ostringstream ss;
			ss << "P:/newpics/big4/zeta";
			ss << imgi;
			//ss << "(";
			//ss << y;
			//ss << ")";
			ss << ".png";
			lodepng::encode(ss.str(), img, w, h);

			imgi += 1;
		}
	}
}












void partialzetasumsanalytic() {
	int offset = 1;
	int width = 1980;
	int height = 1080;
	int xoff = 0;

	int samples = 20;
	int imgi = 0;
	float maxsample = 28.f;
	float scale = 600;
	int granularity = 1;

	std::vector<vec2> traj;
	std::vector<int> traj_meta;
	std::vector<vec2> primal;

	for (int yyy = 0; yyy < 1360000; ++yyy) {
		//int y = nth_prime(yyy);
		//
		int y = yyy;

		float yy = y / float(granularity);

		setwh(width, height);

		gprint(20, 20, "i=%x", yy);

		vec2 prev;
		float total_rotation = 0.f;

		for (int s = 1; s <= samples; ++s) {
			float rotational_contribution = yy * (360.f / s);
			total_rotation += rotational_contribution;

			vec2 offset = vec2(1, 0).rotate(rotational_contribution, vec2()).set_length(nth_length(s));

			//if (!(s % 2)) {
			//	offset = -offset;
			//}

			vec2 newpos = prev + offset;

			//if(outputiter)
			Line(xcoord(prev.x), prev.y*scale + height / 2, xcoord(newpos.x), newpos.y*scale + height / 2, 255);

			prev = newpos;
		}

		auto ffi = get_full_factor_info(y / granularity);
		int num_factors = ffi.prime_factors_sequence.size();

		if (y % granularity == 0 && num_factors == 1) {
			primal.push_back(prev);
		}

		traj.push_back(prev);
		traj_meta.push_back(num_factors);

		prev = traj[0];

		for (int i = 0; i < traj.size(); ++i) {
			auto t = traj[i];
			Line(xcoord(prev.x), prev.y*scale + height / 2, xcoord(t.x), t.y*scale + height / 2, 0, 255, 0);
			//setpix(xcoord(t.x), t.y*scale + height / 2, 0, 255.f, 0);
			prev = t;
		}

		if (primal.size() > 1) {
			prev = primal[0];

			for (size_t i = 1; i < primal.size(); ++i) {
				auto p = primal[i];
				Line(xcoord(prev.x), prev.y*scale + height / 2, xcoord(p.x), p.y*scale + height / 2, 255, 0, 0);
				//setpix(xcoord(p.x), p.y*scale + height / 2, 255, 0, 0);
				prev = p;
			}
		}

		makecross(xcoord(0), height / 2, 0, 0, 255);
		makecross(xcoord(1), height / 2, 0, 0, 255);

		std::ostringstream ss;
		ss << "P:/newpics/big4/zeta";
		ss << imgi;
		//ss << "(";
		//ss << y;
		//ss << ")";
		ss << ".png";
		lodepng::encode(ss.str(), img, w, h);

		imgi += 1;
	}
}

void partialzetasums() {
	partialzetasumsintegers();

}