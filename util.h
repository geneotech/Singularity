#pragma once
#include <fstream>
#define PId 3.1415926535897932384626433832795
#include <iostream>

#include <complex>

typedef std::complex<long double> cc;

static const int g = 7;
static const double p[g + 2] = { 0.99999999999980993, 676.5203681218851,
-1259.1392167224028, 771.32342877765313, -176.61502916214059,
12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
1.5056327351493116e-7 };

static cc cgamma(cc z)
{
	if (real(z)<0.5) {
		return cc(PId) / (sin(cc(PId)*z)*cgamma(cc(1.0) - z));
	}
	z -= 1.0;
	cc x = p[0];
	for (int i = 1; i<g + 2; i++) {
		x += cc(p[i]) / (z + cc(i, 0));
	}
	cc t = z + cc(g + 0.5);
	return cc(sqrt(2 * PId)) * pow(t, z + cc(0.5)) * exp(-t) * x;
}


static cc zeta_alternating(float xx, float yy, const int iters) {
	cc expo(xx, yy);
	cc mult = cc(1, 0) - (cc(2, 0) / pow(cc(2, 0), expo));

	cc rhs;
	for (int i = 1; i < iters; ++i) {
		auto term = cc(1, 0) / pow(cc(i, 0), expo);

		if (i % 2 == 0)
			rhs -= term;
		else
			rhs += term;
	}

	return rhs;// / mult;
}

static cc zeta(float xx, float yy, const int iters = 6260);

static cc zeta(cc s, const int iters = 6260) {
	return zeta(real(s), imag(s), iters);
}

static cc zeta(float xx, float yy, const int iters) {
	cc s(xx, yy);

	//if (xx < 0.5f) {
	//	return std::pow(2, s) * std::pow(PId, s - cc(1)) * std::sin(cc(PId / 2) * s) * cgamma(cc(1) - s) * zeta(cc(1) - s);
	//}
	//else
		return zeta_alternating(xx, yy, iters);
}

static cc spiral(float xx, float yy) {
	cc s(xx, yy);
	return std::pow(2, s) * std::pow(PId, s - cc(1)) * std::sin(cc(PId / 2) * s) * cgamma(cc(1) - s);
}

template<class... Args>
static void Line(float x1, float y1, float x2, float y2, Args... args)
{
	if (!inbounds(x1, y1) && !inbounds(x2, y2)) return;

	// Bresenham's line algorithm
	const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
	if (steep)
	{
		std::swap(x1, y1);
		std::swap(x2, y2);
	}

	if (x1 > x2)
	{
		std::swap(x1, x2);
		std::swap(y1, y2);
	}

	const float dx = x2 - x1;
	const float dy = fabs(y2 - y1);

	float error = dx / 2.0f;
	const int ystep = (y1 < y2) ? 1 : -1;
	int y = (int)y1;

	const int maxX = (int)x2;

	for (int x = (int)x1; x<maxX; x++)
	{
		if (steep)
		{
			setpix(y, x, args...);
		}
		else
		{
			setpix(x, y, args...);
		}

		error -= dy;
		if (error < 0)
		{
			y += ystep;
			error += dx;
		}
	}
}
#include "math/vec2.h"
template<class... Args>
static void makecross(int x, int y, Args... args) {
	Line(x - 10, y, x + 10, y, args...);
	Line(x, y - 10, x, y + 10, args...);
}


#include "typesafe_sprintf.h"
#include "glyphs.h"

template<class... Args>
static void gprint(int x, int y, Args... args) {
	std::string ss = typesafe_sprintf(args...);

	for (auto l : ss) {
		if (l == ' ')
			l = '_';

		auto& g = glyphs[l];

		for (int i = 0; i < g.w; ++i) {
			for (int j = 0; j < g.h; ++j) {
				if (g.img[(j * g.w + i) * 4] == 0)
					setpix(x + i, y + g.h - j - 1, 255);
			}
		}

		x += g.w;
	}
}