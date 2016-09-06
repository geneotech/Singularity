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

#include "../util.h"

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
	
		 int samples = 340;
		 int imgi = 0;
		 float maxsample = 28.f;
		 float scale = 150;
		 
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
		 /*
		 for (float y = from; y < to; y += 0.1f) {
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
		 
		 	for (float xx = 0.5f; xx <= 8.f; xx += 0.05)
		 	{
		 		float yy = y;

		 		auto result = zeta(xx, yy);
		 
		 		int new_x = result.real()*scale + width / 2;
		 		int new_y = result.imag()*scale + height / 2;
		 
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
			
			makecross(critline.real()*scale + width / 2, critline.imag()*scale + height / 2, 255);
			makecross(critbound0.real()*scale + width / 2, critbound0.imag()*scale + height / 2, 255, 0, 255);
			makecross(critbound1.real()*scale + width / 2, critbound1.imag()*scale + height / 2, 255, 0, 255);
			
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
		 }*/
		 // zeta euler spirals

		 for (double y = 0.; y < 450.; y += 0.05)
		 //float y = 60.83177;
		 {

			 auto critline = zeta(0.5f, y);
			 auto critbound0 = zeta(0.0f, y);
			 auto critbound1 = zeta(1.0f, y);

			 //for (float xx = 0.5f; xx <= 2.f; xx += 0.01)
			 float xx = 0.5f;
			 {
				 setwh(width, height);

				 Line(width / 2 - 10, height / 2, width / 2 + 10, height / 2, 255, 0, 0);
				 Line(width / 2, height / 2 - 10, width / 2, height / 2 + 10, 255, 0, 0);

				 {
					 float yy = y;

					 cc prevcomplex;

					 int prev_x = width / 2;
					 int prev_y = height / 2;

					 bool wasprev = true;

					 cc expo(xx, yy);
					 cc mult = cc(1, 0) - (cc(2, 0) / pow(cc(2, 0), expo));

					 cc rhs;

					 for (int i = 1; i <= 600; ++i) {
						 auto term = cc(1, 0) / pow(cc(i, 0), expo);

						 if (i % 2 == 0)
							 rhs -= term;
						 else
							 rhs += term;

						 auto result = rhs;// / mult;

						 int new_x = result.real()*scale + width / 2;
						 int new_y = result.imag()*scale + height / 2;

						 if (wasprev) {
							 //if (c > 0 && (maxsample / samples) * c > 0.5f && (maxsample / samples) * (c - 1) <= 0.5f)
							 //	Line(prev_x, prev_y, new_x, new_y, 255);
							 //else if (c > 0 && (maxsample / samples) * c > 1.0f && (maxsample / samples) * (c - 1) <= 1.0f)
							 //	Line(prev_x, prev_y, new_x, new_y, 0, 0, 255);
							 //else if (c == 0)
							 //	Line(prev_x, prev_y, new_x, new_y, 0, 0, 255);
							 //else
							 //Line(prev_x, prev_y, new_x, new_y, 0, 255, 0);
						 }
						 wasprev = true;
						 //std::cout << (result - prevcomplex) << std::endl;
						 prevcomplex = result;

						 prev_x = new_x;
						 prev_y = new_y;
					 }
				 }
				 
				 {
					 long double param1 = y;
					 typedef vec2t<long double> vd;

					 vd pos(1, 0), prevpos;
					 
					 Line(prevpos.x*scale + width / 2, prevpos.y*scale + height / 2, pos.x*scale + width / 2, pos.y*scale + height / 2, 255);
					 
					 long double dt = 1;

					 for (long double i = 1 + dt; i <= 600; i += dt) {
						 vd lastvec = prevpos - pos;
						 
						 augs::rotate_rad(lastvec, -param1 * (log(i) - log(i - dt)));
						 
						 lastvec.set_length(1.0 / sqrt(i));
						 
						 prevpos = pos;
						 pos += lastvec;

						 Line(prevpos.x*scale + width / 2, prevpos.y*scale + height / 2, pos.x*scale + width / 2, pos.y*scale + height / 2, 0, 255,0);
					 }
				 }
				 /*
				 {
					 long double param1 = y*1000;
					 typedef vec2t<long double> vd;
					 typedef long double ld;

					 vd pos(1, 0), prevpos;
					 vd vel(1, 0);
					 ld angvel = param1;
					 
					 //augs::rotate_rad(vel, -param1);

					 long double dt = 0.1;

					 for (long double i = 1+dt; i <= 600; i += dt) {
						 // integrate

						 //angvel += 0.00001;
						 angvel = -param1 * (log(i) - log(i-dt)) + PId*dt;

						 //vel.set(cos(rads), sin(rads)) * 1.0/sqrt(i);
						 vel.set_length(1.0 / sqrt(i));
						 augs::rotate_rad(vel, angvel * dt);

						 prevpos = pos;
						 pos += vel * dt;

						 Line(prevpos.x*scale + width / 2, prevpos.y*scale + height / 2, pos.x*scale + width / 2, pos.y*scale + height / 2, 255);
					 }
				 }*/


				 {
					 typedef vec2t<long double> vd;
					 typedef long double ld;

					 vd pos(1, 0), prevpos;
					 vd vel(1, 0);
					 ld angvel = PId/2;
					 //scale = 40;

					 //augs::rotate_rad(vel, -param1);

					 long double dt = 0.01;

					 for (long double i = 1 + dt; i <= 600; i += dt) {
						 // integrate
						 
						 //angvel -= 0.04;
						 //vel.set(cos(rads), sin(rads)) * 1.0/sqrt(i);
						 //vel.set_length(1.0 / sqrt(i));
						 augs::rotate_rad(vel, angvel * dt);
						 vel.set_length( 1/sqrt(i));
						 if (i > 300)
							 angvel -= 0.01;
							 //vel.set_length(8 / sqrt(i));

						 prevpos = pos;
						 pos += vel * dt;

						 Line(prevpos.x*scale + width / 2, prevpos.y*scale + height / 2, pos.x*scale + width / 2, pos.y*scale + height / 2, 255);
					 }
				 }

				 //makecross(critline.real()*scale + width / 2, critline.imag()*scale + height / 2, 255);
				 //makecross(critbound0.real()*scale + width / 2, critbound0.imag()*scale + height / 2, 255, 0, 255);
				 //makecross(critbound1.real()*scale + width / 2, critbound1.imag()*scale + height / 2, 255, 0, 255);

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