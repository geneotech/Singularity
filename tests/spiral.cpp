#include "tests.h"
#include "../img.h"

struct vec2 {
	int x;
	int y;
};

int spiral_pos_to_index(vec2 p)
{
	float x = abs(p.x);
	float y = abs(p.y);
	bool q = x > y;

	x = q ? x : y;
	y = q ? p.x + p.y : p.x - p.y;
	y = abs(y) + 4. * x * x + 1.;
	x *= 2.;

	return q ? (p.x > 0. ? y - x - x : y)
		: (p.y > 0. ? y - x : y + x);
}

int ulam_get_map(int x, int y, int n)
{
	x -= (n - 1) / 2;
	y -= n / 2;

	int mx = abs(x), my = abs(y);
	int l = 2 * std::max(mx, my);
	int d = y > x ? l * 3 + x + y : l - x - y;

	return pow(l - 1, 2) + d;
}

#include <cassert>
#include <iostream>
using namespace std;
vec2 spiral_index_to_pos(int i) {
	double sidelen = floor(sqrt(i));
	
	//std::cout << i << " " << sidelen << std::endl;

	double layer = ((sidelen + 1) / 2);
	double offset = i - sidelen*sidelen;
	double segment = (offset / (sidelen + 1));
	double offset2 = offset - 4 * segment + 1 - layer;

	int x, y;
	switch (int(segment)) {
	case 0: x = layer; y = offset2; break;
	case 1: x = -offset2; y = layer; break;
	case 2: x = -layer; y = -offset2; break;
	case 3: x = offset2; y = -layer; break;
	default: assert(0); break;
	}

	return{ x, y };
}
#include "../primes.h"

#include "../modulity.h"

#define SIDE 10000
#define COUNT SIDE*SIDE

int linearized_divisor_plot[COUNT];
int special_col[COUNT];

void arithm_spiral() {
	setwh(SIDE, SIDE);

	int dir = 0;
	int x= SIDE /2, y= SIDE /2;

	for (int n = 1; n < 9500; ++n) {
		bool ispri = is_prime_lookup[n];

		for (int k = 1; k <= n; ++k) {
			auto val = gcd_length(n, k);
			
			if(val)
				setpix(x, y, 255/val);

			// if(val && ispri)
			// 	setpix(x, y, 255, 0, 0);

			if (dir == 0) ++x;
			if (dir == 1) ++y;
			if (dir == 2) --x;
			if (dir == 3) --y;
		}
		++dir;
		dir %= 4;
	}

	lodepng::encode("pics/arithm_spiral.png", img, w, h);
}

void spiral() {
	std::fill(linearized_divisor_plot, linearized_divisor_plot + COUNT, 0);
	std::fill(special_col, special_col + COUNT, 0);
	
	setwh(SIDE, SIDE);
	int n = SIDE;
	
	int l = 0;
	for (int x = 0; l < COUNT; ++x) {
		bool ispri = is_prime_lookup[x];
		//if (!ispri)
		//	continue;
		//	special_col[l] = true;

		if (x % 2 == 0) {
			for (int y = 2; y < x && l < COUNT; ++y) {
				linearized_divisor_plot[l] = gcd_length(x, y) - 1;
				// special_col[l] = ispri;
				++l;
			}
		}
		else {
			for (int y = x-1; y > 1 && l < COUNT; --y) {
				linearized_divisor_plot[l] = gcd_length(x, y) - 1;
				// special_col[l] = ispri;
				++l;
			}
		}
	}

	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			int z = ulam_get_map(y, x, n);

			//if (linearized_divisor_plot[z])
			//	setpix(x, y, 255);
			auto val = linearized_divisor_plot;
			if (val[z]) {
				if (special_col[z] == 0) {
					setpix(x, y, 255 / val[z]);
				}
				//if (special_col[z] == 1){
				//	setpix(x, y, 255, 0, 0);
				//}
			}
		}
	}

	lodepng::encode("pics/ulam.png", img, w, h);
}