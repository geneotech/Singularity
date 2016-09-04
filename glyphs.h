#pragma once

struct glyph {
	std::vector<unsigned char> img;
	unsigned w, h;
};

extern glyph glyphs[256];