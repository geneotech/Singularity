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
	int width = 1920;
	int height = 1080;
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

	for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
		for (int i = 1; i <= std::min(n, height); ++i) {
			double x = n;
			double y = i;
			x /= 100;
			y /= 100;

			//double val = 0.0;

			//{
			//	auto sinval = sin(3.1415926535897932384626433832795 * ((x/y) + x / ((x / (x / y + 1)) + 1) + x / ((x/y)+1) ));
			//	val += 1 - sinval*sinval;
			//}

			//{
			//	auto sinval = sin(3.1415926535897932384626433832795 * x / y);
			//	val += 1 - sinval*sinval;
			//}
			//
			//if (val > 0)
			//	setpix(j - offset, i, 255 * val);


			//auto val = modulity(n, i) + 1;
			auto a = sin(PId*x);
			auto b = sin(PId * y / x);
			auto val = a*a + b*b;
			//val += 2;

			if (val != 0)
				setpix(j - offset, i, val*100);
		}
	}
	lodepng::encode("newpics/sympho.png", img, w, h);

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