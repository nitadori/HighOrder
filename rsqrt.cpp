// To compile in MacOS:
// g++ -O2 -Wall -mfma -Wa,-q --std=c++11 rsqrt.cpp -I . libqd.a

#include <iostream>
#include <cstdlib>
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
#include <qd/qd_real.h>
#include <qd/qd_inline.h>

dd_real rsqrt1(const dd_real &x){
	double y0 = 1.0 / sqrt(to_double(x));
	return dd_real(y0);
}

dd_real rsqrt2(const dd_real &x){
	//  d  * d  -> dd
	//  dd * dd -> dd
	//  d  - dd -> d
	//  d  * dd -> dd
	//  d  + d  -> dd
	double y0 = 1.0 / sqrt(to_double(x));
	double h  = (1.0 - x * dd_real::sqr(y0)).x[0];
	// dd_real y1 = y0 + y0 * mul_pwr2(h, 0.5);
	dd_real y1 = dd_real::add(y0,  (0.5 * y0) * h);
	return y1;
}

dd_real rsqrt3(const dd_real &x){
	//  d  * d  -> dd
	//  dd * dd -> dd
	//  d  - dd -> dd
	//  d  * dd -> dd
	double y0 = 1.0 / sqrt(to_double(x));
	dd_real c = 3.0 - x * dd_real::sqr(y0);
	dd_real y1 = (0.5 * y0) * c;
	return y1;
}

dd_real rsqrt4(const dd_real &x){
	//  d  * d  -> dd
	//  dd * dd -> dd
	//  d  - dd -> dd
	//  d  + dd -> dd
	//  d  * dd -> dd
	double y0 = 1.0 / sqrt(to_double(x));
	dd_real h = 1.0 - x * dd_real::sqr(y0);
	dd_real c = dd_real(0.5, (3./8.) * to_double(h));
	// dd_real y1 = y0 + y0 * h * c;
	dd_real y1 = y0 * (1.0 +  h * c);
	return y1;
}

dd_real rsqrt5(const dd_real &x){
	double y0 = 1.0 / sqrt(to_double(x));
	double h = to_double(1.0 - x * dd_real::sqr(y0));
	dd_real y1 = dd_real::add( y0,  y0 * (h * (0.5 + h * (3./8.))) );
	return y1;
}

dd_real rsqrt_qd(const dd_real &x){
	return to_dd_real(1.0 / sqrt(qd_real(x)));
}


template <typename F>
void accruracy_check(const dd_real &x, F fun){
	dd_real y = fun(x);
	dd_real h = 1.0 - x * sqr(y);

	qd_real yq = 1.0 / sqrt(qd_real(x));
	double qerr = to_double(y - yq) * sqrt(to_double(x));
	double derr = to_double(y - to_dd_real(yq)) * sqrt(to_double(x));

	printf("%A %A\n", y.x[0], y.x[1]);

	std::cout << y << " " << h << " " << qerr << " " << derr << std::endl;
}

int main(){
	using std::cout;
	using std::endl;

	srand(20170728);

	for(int i=0; i<10; i++){

		dd_real x = ddrand();

		cout.precision(32);

		cout << x << endl;
		cout << std::scientific << to_double(x) << endl;

		/*
		dd_real y = 1.0 / sqrt(x);
		dd_real h = 1.0 - x * sqr(y);

		cout << y << " " << h << endl;
		*/
		// accruracy_check(x, rsqrt1); 
		accruracy_check(x, [](const dd_real &x){ return 1.0 / sqrt(x); }); 
		accruracy_check(x, rsqrt2); 
		// accruracy_check(x, rsqrt3); 
		// accruracy_check(x, rsqrt4); 
		accruracy_check(x, rsqrt5); 
		accruracy_check(x, rsqrt_qd); 
	}
}
