#pragma once
#include <vector>
#include <algorithm>
#include "lodepng.h"
#include <sstream>

static int w = 3302;
static int h = 400;
static std::vector<unsigned char> img;

static void setpix(int x, int y, int val) {
	if (x >= w || y >= h || y < 0) return;

	img[(w * (h - y - 1) + x) * 4] = val;
	img[(w * (h - y - 1) + x) * 4 + 1] = val;
	img[(w * (h - y - 1) + x) * 4 + 2] = val;
}

static void setpix(int x, int y, int r, int g, int b) {
	if (x >= w || y >= h) return;
	img[(w * (h - y - 1) + x) * 4] = r;
	img[(w * (h - y - 1) + x) * 4 + 1] = g;
	img[(w * (h - y - 1) + x) * 4 + 2] = b;
}

static void setpixrgb(int x, int y, unsigned RGBint) {
	if (x >= w || y >= h) return;
	img[(w * (h - y - 1) + x) * 4] = RGBint & 255;
	img[(w * (h - y - 1) + x) * 4 + 1] = (RGBint >> 8) & 255;
	img[(w * (h - y - 1) + x) * 4 + 2] = (RGBint >> 16) & 255;
}

static void setwh(int _w, int _h) {
	w = _w;
	h = _h;

	img.clear();
	img.resize(w * h * 4, 0);

	for (int i = 0; i < w*h; ++i) {
		img[i * 4 + 3] = 255;
	}
}

#define COLS 7
static void setcol(int x, int y, int c) {
	switch (c) {
	case 0: setpix(x, y, 255); break;
	case 1: setpix(x, y, 255, 0, 0); break;
	case 2: setpix(x, y, 255, 255, 0); break;
	case 3: setpix(x, y, 0, 255, 0); break;
	case 4: setpix(x, y, 0, 0, 255); break;
	case 5: setpix(x, y, 0, 255, 255); break;
	case 6: setpix(x, y, 255, 0, 255); break;
	}
}