#include "tests/tests.h"
#include "img.h"
#include "modulity.h"

#include "math/vec2.h"

#define PId 3.1415926535897932384626433832795

double analytic_floor(double x) {
	double fracpart = 0.0;
	fracpart = 1 / 2.0;

	for (int m = 1; m < 30; ++m) {
		fracpart += ((-1) / (m*PId)) * sin(2*PId*m*x);
	}

	return x - fracpart;
}

double eps = 0.0001;
#include <iostream>

void hyperbolalattice() {
	int offset = 1;
	int width = 1980;
	int height = 1080;

	double pixels_per_square = 100;

	double interval = 0.0002;

	setwh(width, height);

	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
		}
	}

	//for (double x = interval; x < width; x += interval) {
	//	double val = analytic_floor(x);
	//	setpix(x* pixels_per_square, val* pixels_per_square, 255);
	//}
	//std::ostringstream ss;
	//ss << "pics/analytic";
	//ss << ".png";
	//lodepng::encode(ss.str(), img, w, h);

	double dd[4];
	vec2 vv[4];

	double from = 2;
	double to = 13;

	int maxs = 3;

	for (double m = from; m <= to; m += 0.05) {
		//for (double m = 2; m < 11;m+=0.05) {
		setwh(width, height);
		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < height; ++y) {
				setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
			}
		}

		for (double x = interval; x < width/pixels_per_square; x += interval) {
			const double val = m / x;

			vec2 p(x, val);
			vec2 ip(floor(p.x), floor(p.y));

			double dist = 0.0;
			
			{
				double d = (p - ip).length_sq();

				if (d > eps)
					dist += 1 / d;
			}

			int s = 3;
			bool was_bigger_than_epsilon = false;

			while(true) {
				was_bigger_than_epsilon = false;

				int corner_center_offset = (s - 1) / 2;
				int cco = corner_center_offset;

				for (int ss = 0; ss < s-1; ++ss) {
					vec2 lt(-cco, cco - ss);
					vec2 lb(-cco + ss, -cco);
					vec2 rb(cco, -cco + ss);
					vec2 rt(cco - ss, cco);

					vec2 points[4] = { ip+lt, ip + lb, ip + rb, ip + rt };

					for (auto& pp : points) {
						double d = (pp - p).length_sq();
						d = 1 / d;

						if (d > eps) {
							was_bigger_than_epsilon = true;
						}

						dist += d;
					}
				}

				if (!was_bigger_than_epsilon)
					break;

				s += 2;
			}

			if (s > maxs) {
				maxs = s;
			}

			dist /= 20.0;

			setpix(x* pixels_per_square, dist* pixels_per_square, 255);
			setpix(x* pixels_per_square, val* pixels_per_square, 0, 255, 0);
		}
		std::ostringstream ss;
		ss << "pics/latticeseries";
		ss << std::fixed << m;
		ss << ".png";
		lodepng::encode(ss.str(), img, w, h);
	}


	std::cout << "Maximum side: " << maxs;
	int abca;
	std::cin >> abca;

	return;

	for (double m = from; m <= to; m += 0.05) {
		//for (double m = 2; m < 11;m+=0.05) {
		setwh(width, height);
		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < height; ++y) {
				setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
			}
		}
	
		for (double x = interval; x < width; x += interval) {
			const double val = m / x;
	
			vec2 lb(floor(x), floor(val));
			vec2 lt(lb.x, lb.y + 1);
			vec2 rb(lb.x + 1, lb.y);
			vec2 rt(lb.x + 1, lb.y + 1);
	
			vec2 p(x, val);
	
			vv[0] = lb;
			vv[1] = lt;
			vv[2] = rb;
			vv[3] = rt;
			dd[0] = (p - lb).length();
			dd[1] = (p - lt).length();
			dd[2] = (p - rb).length();
			dd[3] = (p - rt).length();

			int mindist = std::min_element(dd + 0, dd + 4) - dd;

			//double dist = dd[0] * dd[1] * dd[2] * dd[3];
			////dist *= (dd[0] + dd[1] + dd[2] + dd[3]);
			//dist *= 20;


			vec2 c = lb;// vv[mindist];

			vec2 q[8] = {
				c + vec2(1, 0),
				c + vec2(1, -1),
				c + vec2(0, -1),
				c + vec2(-1, -1),
				c + vec2(-1, 0),
				c + vec2(-1, 1),
				c + vec2(0, 1),
				c + vec2(1, 1)
			};

			double product = (p - c).length();

			for (auto& qq : q) {
				product *= (p - qq).length();
			}

			dd[0] = std::max(dd[0], eps);
			dd[1] = std::max(dd[1], eps);
			dd[2] = std::max(dd[2], eps);
			dd[3] = std::max(dd[3], eps);

			double dist = product;
			//dist = dd[0] * dd[1] * dd[2] * dd[3] * 10;
			dist = (1/dd[0] + 1 / dd[1] + 1 / dd[2] + 1 / dd[3]);
			setpix(x* pixels_per_square, dist* pixels_per_square, 255);
			setpix(x* pixels_per_square, val* pixels_per_square, 0, 255, 0);
		}
		std::ostringstream ss;
		ss << "pics/latticefloor";
		ss << std::fixed << m;
		ss << ".png";
		lodepng::encode(ss.str(), img, w, h);
	}
}