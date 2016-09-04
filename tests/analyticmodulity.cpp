#include "tests.h"
#include "../img.h"
#include "../modulity.h"

long long
_MM(long long a, long long b)
{
	while (true) {
		if (!(a%b)) {
			return b;
		}
		b = a%b;
	}
}

#include <fstream>
#define PId 3.1415926535897932384626433832795
#include <iostream>

#include "../primes.h"
#include "../factorizator.h"
#include <complex>

typedef std::complex<long double> cc;

static const int g = 7;
static const double p[g + 2] = { 0.99999999999980993, 676.5203681218851,
-1259.1392167224028, 771.32342877765313, -176.61502916214059,
12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
1.5056327351493116e-7 };

cc cgamma(cc z)
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


cc zeta_alternating(float xx, float yy, const int iters) {
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

	return rhs / mult;
}

cc zeta(float xx, float yy, const int iters = 6260);

cc zeta(cc s, const int iters = 6260) {
	return zeta(real(s), imag(s), iters);
}

cc zeta(float xx, float yy, const int iters) {
	cc s(xx, yy);

	if (xx < 0.5f) {
		return std::pow(2, s) * std::pow(PId, s - cc(1)) * std::sin(cc(PId / 2) * s) * cgamma(cc(1) - s) * zeta(cc(1) - s);
	}
	else
		return zeta_alternating(xx, yy, iters);
}

template<class... Args>
void Line( float x1,  float y1,  float x2,  float y2, Args... args)
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
#include "../math/vec2.h"
template<class... Args>
void makecross(int x, int y, Args... args) {
	Line(x - 10, y, x + 10, y, args...);
	Line(x, y- 10, x, y + 10, args...);
}


#include "../typesafe_sprintf.h"
#include "../glyphs.h"

template<class... Args>
void gprint(int x, int y, Args... args) {
	std::string ss = typesafe_sprintf(args...);

	for (auto l : ss) {
		if (l == ' ')
			l = '_';

		auto& g = glyphs[l];

		for (int i = 0; i < g.w; ++i) {
			for (int j = 0; j < g.h; ++j) {
				if(g.img[(j * g.w + i) * 4] == 0)
					setpix(x + i, y + g.h - j - 1, 255);
			}
		}

		x += g.w;
	}
}

void analyticmodulity() {

	//for (int i = 2; i < 100; ++i) {
	//	cout << i << ": ";
	//	
	//	int p = i;
	//	while (true) {
	//		if (is_prime_lookup[p])
	//			break;
	//		
	//		--p;
	//	}
	//	
	//	cout << "M(" << i << ", " << p << ") = " << modulity(i, p) << endl;
	//}
	//
	//std::cout << "Done.";
	//int ab;
	//cin >> ab;


	int offset = 1;
	int width = 1280;
	int height = 720;
	int xscale = 1;

	//_MM(13, 70);
	//gcd_length(13, 70);
	//std::cout << 27 % 8;
	////std::cin >> xscale;
	//
	setwh(width, height);
	//
	//{
	//	std::ofstream gen_file("gcds.txt");
	//
	//	for (int i = 1; i < 1000; ++i) {
	//		for (int j = 1; j < 1000; ++j) {
	//			if(i >= j && gcd_length(i, j) != _MM(i, j)+1)
	//				gen_file << i << " " << j << " " << gcd_length(i, j) << " " << _MM(i, j) << std::endl;
	//		}
	//	}
	//}

	typedef std::complex<long double> cc;

	// zeta three comparatory spirals

		//float from = 200.0;
		//float to = 210.f;
		//int imgi = 0;
		//for (float y = from; y < to; y += 0.01f) {
		//
		//	std::ostringstream fn;
		//	fn << "testlens/lengths";
		//	fn << imgi;
		//	fn << "(";
		//	fn << y;
		//	fn << ")";
		//	fn << ".txt";
		//
		//	std::ofstream gen_file(fn.str());
		//	//float xx = y;
		//	//float yy = 14.134725141734f;
		//
		//	setwh(width, height);
		//
		//	for(int c = 0; c < 3; ++c)
		//	{
		//
		//		float xx = 0.f;
		//
		//		if (c == 0)
		//			xx = 0.25f;
		//		if (c == 1)
		//			xx = 0.5f;
		//		if (c == 2)
		//			xx = 0.75f;
		//
		//		float yy = y;
		//
		//		cc expo(xx, yy);
		//		auto mult = cc(1, 0) - (cc(2, 0) / pow(cc(2, 0), expo));
		//
		//		int iters = 160;
		//
		//	int prev_x = width / 2;
		//	int prev_y = height / 2;
		//	
		//	cc rhs;
		//
		//	cc prevcomplex;
		//	cc pprevcomplex;
		//	for (int i = 1; i < iters; ++i) {
		//		auto term = cc(1, 0) / pow(cc(i, 0), expo);
		//
		//		if (i % 2 == 0)
		//			rhs -= term;
		//		else
		//			rhs += term;
		//
		//		auto result = rhs / mult;
		//
		//		int new_x = result.real()*150.f + width / 2;
		//		int new_y = result.imag()*150.f + height / 2;
		//
		//		if(c == 0)
		//			Line(prev_x, prev_y, new_x, new_y, 255);
		//		if (c == 1)
		//			Line(prev_x, prev_y, new_x, new_y, 255, 0, 0);
		//		if (c == 2)
		//			Line(prev_x, prev_y, new_x, new_y, 0, 255, 0);
		//
		//		vec2 ppr(pprevcomplex.real(), pprevcomplex.imag());
		//		vec2 f1(prevcomplex.real(), prevcomplex.imag());
		//		vec2 f2(result.real(), result.imag());
		//
		//		//gen_file << (f2 - f1).length() << std::endl;
		//		//gen_file << (f2 - f1).radians_between(ppr - f1) << std::endl;
		//		gen_file << result.real() << " " << result.imag() << std::endl;
		//		pprevcomplex = prevcomplex;
		//		prevcomplex = result;
		//
		//		prev_x = new_x;
		//		prev_y = new_y;
		//	}
		//	}
		//
		//	std::ostringstream ss;
		//	ss << "newpics2/zeta";
		//	ss << imgi;
		//	ss << "(";
		//	ss << y;
		//	ss << ")";
		//	ss << ".png";
		//	lodepng::encode(ss.str(), img, w, h);
		//
		//	imgi += 1;
		//}

		// zeta trajectory along strip
		// single image samples all critical x for given y
	
		 float from = -1;
		 float to = 100;
		 int samples = 240;
		 int imgi = 0;
		 float maxsample = 8.f;
		 
		 std::cout << "zeta(-7.f, 7.f) = " << zeta(-7.f, 7.f) << std::endl;
		 std::cout << "zeta(1, 7.f)    = " << zeta(1, 7.f) << std::endl;
		 std::cout << "zeta(0, 7.f)    = " << zeta(0, 7.f) << std::endl;
		 std::cout << "zeta(0.1, 7.f)  = " << zeta(0.1, 7.f) << std::endl;
		 std::cout << "zeta(0.2, 7.f)  = " << zeta(0.2, 7.f) << std::endl;
		 std::cout << "zeta(0.3, 7.f)  = " << zeta(0.3, 7.f) << std::endl;
		 std::cout << "zeta(0.4, 7.f)  = " << zeta(0.4, 7.f) << std::endl;
		 std::cout << "zeta(0.5, 7.f)  = " << zeta(0.5, 7.f) << std::endl;
		 std::cout << "zeta(-1, 7.f)   = " << zeta(-1, 7.f) << std::endl;
		 std::cout << "zeta(-2, 7.f)   = " << zeta(-2, 7.f)  << std::endl;

		 for (float y = from; y < to; y += 0.05f) {
		 	setwh(width, height);
		 	
		 	Line(width / 2 - 10, height / 2 , width / 2 + 10, height / 2 , 255, 0, 0);
		 	Line(width/2 , height/2 -10, width / 2 , height / 2 + 10, 255, 0, 0);
		 
		 	cc prevcomplex;
		 
		 	int prev_x = width / 2;
		 	int prev_y = height / 2;
		 
		 	bool wasprev = false;

			auto critline = zeta(0.5f, y);
			auto critbound0 = zeta(0.0f, y);
			auto critbound1 = zeta(1.0f, y);
		 
		 	for (int c = -samples*4; c < samples; ++c)
		 	{
		 		float xx = (maxsample /samples) * c;
		 		float yy = y;
		 
		 		auto result = zeta(xx, yy);
		 
		 		int new_x = result.real()*150.f + width / 2;
		 		int new_y = result.imag()*150.f + height / 2;
		 
		 		if (wasprev) {
		 			//if (c > 0 && (maxsample / samples) * c > 0.5f && (maxsample / samples) * (c - 1) <= 0.5f)
		 			//	Line(prev_x, prev_y, new_x, new_y, 255);
					//else if (c > 0 && (maxsample / samples) * c > 1.0f && (maxsample / samples) * (c - 1) <= 1.0f)
					//	Line(prev_x, prev_y, new_x, new_y, 0, 0, 255);
					//else if (c == 0)
					//	Line(prev_x, prev_y, new_x, new_y, 0, 0, 255);
		 			//else
		 				Line(prev_x, prev_y, new_x, new_y, 0, 255, 0);
		 		}
		 		wasprev = true;
		 
		 		prevcomplex = result;
		 
		 		prev_x = new_x;
		 		prev_y = new_y;
		 	}
			
			makecross(critline.real()*150.f + width / 2, critline.imag()*150.f + height / 2, 255);
			makecross(critbound0.real()*150.f + width / 2, critbound0.imag()*150.f + height / 2, 255, 0, 255);
			makecross(critbound1.real()*150.f + width / 2, critbound1.imag()*150.f + height / 2, 255, 0, 255);
			
			gprint(20, 20, "i=%x", y);

		 	std::ostringstream ss;
		 	ss << "P:/newpics/big2/zeta";
		 	ss << imgi;
		 	//ss << "(";
		 	//ss << y;
		 	//ss << ")";
		 	ss << ".png";
		 	lodepng::encode(ss.str(), img, w, h);
		 
		 	imgi += 1;
		 }
		 
		//zeta trajectory along strip
		//single image samples given x for several y
/*
		 float from = 0;
		 float to = 100.f;
		 int samples = 100;
		 int imgi = 100;
		 float maxsample = 7.f;
		 
		 for (int c = 0; c < samples; ++c) {
		 	setwh(width, height);
		 	
		 	Line(width / 2 - 10, height / 2 , width / 2 + 10, height / 2 , 255, 0, 0);
		 	Line(width/2 , height/2 -10, width / 2 , height / 2 + 10, 255, 0, 0);
		 
		 	cc prevcomplex;
		 
		 	int prev_x = width / 2;
		 	int prev_y = height / 2;
		
			float xx = (maxsample / samples) * c + 1.f;
		 
		 	bool wasprev = false;
			
			for (float y = from; y < to; y += 0.05f)
		 	{
		 		float yy = y;
		 
		 		cc expo(xx, yy);
		 		auto mult = cc(1, 0) - (cc(2, 0) / pow(cc(2, 0), expo));
		 
		 		int iters = 6260;
		 
		 		cc rhs;
		 		for (int i = 1; i < iters; ++i) {
		 			auto term = cc(1, 0) / pow(cc(i, 0), expo);
		 
		 			if (i % 2 == 0)
		 				rhs -= term;
		 			else
		 				rhs += term;
		 		}
		 
		 		auto result = rhs / mult;
		 
		 		int new_x = result.real()*150.f + width / 2;
		 		int new_y = result.imag()*150.f + height / 2;
		 
		 		if (wasprev) {
		 			Line(prev_x, prev_y, new_x, new_y, 255);
		 		}
		 		wasprev = true;
		 
		 		prevcomplex = result;
		 
		 		prev_x = new_x;
		 		prev_y = new_y;
		 	}
		 
		 	std::ostringstream ss;
		 	ss << "P:/newpics/traj/zeta";
		 	ss << imgi;
		 	ss << "(";
		 	ss << xx;
		 	ss << ")";
		 	ss << ".png";
		 	lodepng::encode(ss.str(), img, w, h);
		 
		 	imgi += 1;
		 }
		 */

	//	for (unsigned long long x = 1; x < width; ++x) {
	//		for (unsigned long long y = 1; y < height; ++y) {
	//			int primes = 0;
	//
	//			for (int n = 1; n < 100; ++n) {
	//				if ((y) / n == 0)
	//					break;
	//
	//				int evaluated = 2 * (n + (n+x) / ((y) / n)) - 1;
	//				
	//				if (evaluated < max_primes() && is_prime_lookup[evaluated])
	//					++primes;
	//				//else break;
	//			}
	//
	//			//auto info = get_full_factor_info(val);
	//
	//			auto col = primes * 10;
	//			col = std::min(255, col);
	//
	//			setpix(x, y, col);
	//		}
	//	}
	
	//for (int m = 0; m <= 1/ 0.00005; ++m) {
	//	float expo = 1 + m*0.00005;
	//
	//	for (unsigned long long x = 1; x < width; ++x) {
	//		for (unsigned long long y = 1; y < height; ++y) {
	//			auto val = pow(x, expo) + pow(y, expo);
	//
	//			auto info = get_full_factor_info(val);
	//
	//			//if (val < max_primes() && is_prime_lookup[val])
	//			setpix(x, y, 255 / (info.prime_factors_sequence.size()));
	//		}
	//	}
	//	std::ostringstream ss;
	//	ss << "newpics/sympho";
	//	ss << std::fixed << expo;
	//	ss << ".png";
	//	lodepng::encode(ss.str(), img, w, h);
	//}

	
	//for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
	//	for (int i = 1; i <= height; ++i) {
	//		int ix = n;
	//		int iy = i;
	//		double x = ix;
	//		double y = iy;
	//
	//		//double val = 0.0;
	//
	//		//{
	//		//	auto sinval = sin(3.1415926535897932384626433832795 * ((x/y) + x / ((x / (x / y + 1)) + 1) + x / ((x/y)+1) ));
	//		//	val += 1 - sinval*sinval;
	//		//}
	//
	//		//{
	//		//	auto sinval = sin(3.1415926535897932384626433832795 * x / y);
	//		//	val += 1 - sinval*sinval;
	//		//}
	//		//
	//		//if (val > 0)
	//		//	setpix(j - offset, i, 255 * val);
	//
	//
	//		//auto val = modulity(n, i) + 1;
	//		auto val = (ix*ix*ix + iy*iy*iy);
	//		//val += 2;
	//		if(val > 0 && val < max_primes() && is_prime_lookup[val])
	//			setpix(j - offset, i, 255);
	//
	//		//if (val != 0)
	//		//	setpix(j - offset, i, val);
	//	}
	//}
	//lodepng::encode("newpics/sympho.png", img, w, h);

	//setwh(width, height);
	//
	//for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
	//	for (int i = 1; i <= std::min(n, height); ++i) {
	//		double x = n;
	//		double y = i;
	//
	//		auto val = _MM(n, i);
	//
	//		if (val > 0)
	//			setpix(j - offset, i, 255 / val);
	//	}
	//}
	//lodepng::encode("pics/modulo.png", img, w, h);
}