#include "tests.h"
#include "../img.h"

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

	setwh(width, height);

	for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
		for (int i = 1; i <= std::min(n, height); ++i) {
			double x = n;
			double y = i;

			double val = 0.0;

			//{
			//	auto sinval = sin(3.1415926535897932384626433832795 * ((x/y) + x / ((x / (x / y + 1)) + 1) + x / ((x/y)+1) ));
			//	val += 1 - sinval*sinval;
			//}

			{
				auto sinval = sin(3.1415926535897932384626433832795 * x / y);
				val += 1 - sinval*sinval;
			}

			if (val > 0)
				setpix(j - offset, i, 255 * val);
		}
	}
	lodepng::encode("pics/sympho.png", img, w, h);

	setwh(width, height);

	for (int j = offset, n = 0; j <= width + offset; ++j, n = (j)) {
		for (int i = 1; i <= std::min(n, height); ++i) {
			double x = n;
			double y = i;

			auto val = (n % i) + 1;

			if (val > 0)
				setpix(j - offset, i, 255 / val);
		}
	}
	lodepng::encode("pics/modulo.png", img, w, h);
}