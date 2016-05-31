#include "tests/tests.h"
#include "img.h"
#include "modulity.h"

#include "math/vec2.h"
#include "modulity.h"

#define PId 3.1415926535897932384626433832795

double analytic_floor(double x) {
	double fracpart = 0.0;
	fracpart = 1 / 2.0;

	for (int m = 1; m < 30; ++m) {
		fracpart += ((-1) / (m*PId)) * sin(2*PId*m*x);
	}

	return x - fracpart;
}

double eps = 0.00001;
#include <iostream>


struct plotter {
	bool plotted = false;
	int last_y = -1;

	template<class... Args>
	void plot(int x, int y, Args... args) {
		if (plotted) {
			if (last_y < y - 1) {
				int correction = std::min(18000, y - 1 - last_y);
				for (int i = 1; i <= correction; ++i) {
					setpix(x-1, last_y+i, args...);
				}
			}
			else if (last_y > y + 1) {
				int correction = std::min(18000, last_y - (y+1));
				for (int i = 1; i <= correction; ++i) {
					setpix(x - 1, last_y - i, args...);
				}
			}
		}
		
		setpix(x, y, args...);

		last_y = y;
		plotted = true;
	}
};

void hyperbolalattice() {
	int offset = 1;
	double yoff = 0;
	int width = 1980;
	int height = 1080;

	double pixels_per_square = 90;

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

	double from = 11*17;
	double to = 11 * 17;

	for (double m = from; m <= to; m += 0.05) {
		//for (double m = 2; m < 11;m+=0.05) {
		setwh(width, height);
		
		if (pixels_per_square > 1)
		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < height; ++y) {
				setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
			}
		}

		plotter hyper;
		plotter fmx;

		for (int x = 0; x < width; ++x) {
			double fx = x / pixels_per_square;
			const double val = m / fx;

			vec2 lb(floor(fx), floor(val));
			vec2 lt(lb.x, lb.y + 1);
			vec2 rb(lb.x + 1, lb.y);
			vec2 rt(lb.x + 1, lb.y + 1);

			vec2 p(fx, val);

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
			dist = dd[0] * dd[1] * dd[2] * dd[3] * 10;
			//dist = (1 / dd[0] + 1 / dd[1] + 1 / dd[2] + 1 / dd[3]);

			fmx.plot(x, dist* pixels_per_square, 255);
			hyper.plot(x, val* pixels_per_square, 0, 255, 0);
		}

		std::ostringstream ss;
		ss << "newpics/latticefloor";
		ss << std::fixed << m;
		ss << ".png";
		lodepng::encode(ss.str(), img, w, h);
	}

	int maxs = 3;
	//for (double m = from; m <= to; m += 0.05) {
	//	//for (double m = 2; m < 11;m+=0.05) {
	//	setwh(width, height);
	//	if (pixels_per_square > 1)
	//		for (int x = 0; x < width; ++x) {
	//			for (int y = 0; y < height; ++y) {
	//				setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
	//			}
	//		}
	//
	//	plotter hyper;
	//
	//	for (int x = 0; x < width; ++x) {
	//		const double fx = (x+offset) / pixels_per_square;
	//		const double val = gcdd(fx, 1039/fx);
	//
	//		hyper.plot(x, val* pixels_per_square, 255, 255, 255);
	//	}
	//
	//	std::ostringstream ss;
	//	ss << "newpics/frac";
	//	ss << std::fixed << m;
	//	ss << ".png";
	//	lodepng::encode(ss.str(), img, w, h);
	//}
	//
	//return;

	//for (double m = from; m <= to; m += 0.05) {
	//	//for (double m = 2; m < 11;m+=0.05) {
	//	setwh(width, height);
	//	if(pixels_per_square > 1)
	//	for (int x = 0; x < width; ++x) {
	//		for (int y = 0; y < height; ++y) {
	//			setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
	//		}
	//	}
	//
	//	plotter hyper;
	//
	//	for (int x = 0; x < width; ++x) {
	//		const double fx = x / pixels_per_square;
	//		double val = fx * (1 - sin(fx));
	//
	//		const double fy = m / fx;
	//		const double valy = fy * (1 - sin(PId * fy));
	//
	//		vec2 hyp(fx, m / fx);
	//		vec2 latt(val, valy);
	//		auto xs = sin(PId * m / fx);
	//		auto ys = sin(PId * fx);
	//		//xs = xs*xs;
	//		//ys = ys*ys;
	//
	//		val = xs*xs + ys*ys;
	//		//val = pow(xs, ys) *pow(ys, xs);
	//		val = val*2000;
	//
	//		hyper.plot(x, val* pixels_per_square, 255, 255, 255);
	//	}
	//
	//	std::ostringstream ss;
	//	ss << "newpics/frac";
	//	ss << std::fixed << m;
	//	ss << ".png";
	//	lodepng::encode(ss.str(), img, w, h);
	//}
	//
	//return;

	//for (double m = from; m <= to; m += 0.05) {
	//	//for (double m = 2; m < 11;m+=0.05) {
	//	setwh(width, height);
	//		if(pixels_per_square > 1)
	//	for (int x = 0; x < width; ++x) {
	//		for (int y = 0; y < height; ++y) {
	//			setpix(x * pixels_per_square, y * pixels_per_square, 255, 0, 0);
	//		}
	//	}
	//
	//	plotter hyper;
	//	plotter fmx;
	//
	//	for (int x = 0; x < width; ++x) {
	//		const double fx = (x+ offset )/pixels_per_square;
	//		const double val = m / fx;
	//
	//		vec2 p(fx, val);
	//
	//		//double dist = 0.0;
	//		//
	//		//for (int t = 1; t <= m; ++t) {
	//		//	double d = sin(PId*t / x);
	//		//
	//		//	if (d < eps) {
	//		//		d = eps;
	//		//	}
	//		//
	//		//	dist += 1/d;
	//		//}
	//		//
	//		//dist /= 5000.0;
	//		vec2 ip(floor(p.x), floor(p.y));
	//		
	//		double dist = 0.0;
	//		
	//		bool was_bigger_than_epsilon = true;
	//
	//		{
	//			double d = (p - ip).length();
	//
	//			if (d < eps) {
	//				dist = std::numeric_limits<double>::max() / 2;
	//				was_bigger_than_epsilon = false;
	//			}
	//			else {
	//				d = d*d*d;
	//				dist = 1 / d;
	//			}
	//		}
	//		
	//		int s = 3;
	//		
	//		while(was_bigger_than_epsilon) {
	//			was_bigger_than_epsilon = false;
	//		
	//			int corner_center_offset = (s - 1) / 2;
	//			int cco = corner_center_offset;
	//		
	//			for (int ss = 0; ss < s-1; ++ss) {
	//				vec2 lt(-cco, cco - ss);
	//				vec2 lb(-cco + ss, -cco);
	//				vec2 rb(cco, -cco + ss);
	//				vec2 rt(cco - ss, cco);
	//		
	//				vec2 points[4] = { ip+lt, ip + lb, ip + rb, ip + rt };
	//		
	//				for (auto& pp : points) {
	//					//if (pp.x < 0 || pp.y < 0)
	//					//	continue;
	//
	//					double d = (pp - p).length();
	//
	//					if (d < eps) {
	//						dist = std::numeric_limits<double>::max()/2;
	//						was_bigger_than_epsilon = false;
	//						break;
	//					}
	//					
	//					d = d*d*d;
	//					d = 1 / d;
	//		
	//					if (d > eps) {
	//						was_bigger_than_epsilon = true;
	//					}
	//		
	//					dist += d;
	//				}
	//			}
	//		
	//			if (!was_bigger_than_epsilon)
	//				break;
	//		
	//			s += 2;
	//		}
	//		
	//		if (s > maxs) {
	//			maxs = s;
	//		}
	//		
	//		dist = (sqrt(sqrt(dist)));
	//
	//		fmx.plot(x, (dist+ yoff)* pixels_per_square, 255);
	//		hyper.plot(x, (val+ yoff)* pixels_per_square, 0, 255, 0);
	//	}
	//	std::ostringstream ss;
	//	ss << "pics/trial/illustration";
	//	ss << std::fixed << m;
	//	ss << ".png";
	//	lodepng::encode(ss.str(), img, w, h);
	//}
	//
	//
	//std::cout << "Maximum side: " << maxs;
	//int abca;
	//std::cin >> abca;
	//
	//return;

}