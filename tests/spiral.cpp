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
#include <cassert>
vec2 spiral_index_to_pos(int i) {
	int sidelen = sqrt(i);
	int layer = floor((sidelen + 1) / 2);
	int offset = i - sidelen*sidelen;
	int segment = floor(offset / (sidelen + 1));
	int offset2 = offset - 4 * segment + 1 - layer;

	int x, y;
	switch (segment) {
	case 0: x = layer; y = offset2; break;
	case 1: x = -offset2; y = layer; break;
	case 2: x = -layer; y = -offset2; break;
	case 3: x = offset2; y = -layer; break;
	default: assert(0); break;
	}

	return{ x, y };
}

void spiral() {

}