#include <iostream>
#include <cstdlib>
#include <cassert>
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
#include <qd/qd_real.h>
#include <qd/qd_inline.h>

dd_real mult_sub1(const dd_real &a, const dd_real &b){
	double p0, p1, p2, p3;
	double q0, q1, q2;

	p0 = qd::two_prod(a.x[0], b.x[0], q0);

	p1 = qd::two_prod(a.x[0], b.x[1], q1);
	p2 = qd::two_prod(a.x[1], b.x[0], q2);

	p3 = a.x[1] * b.x[1];


	p0 -= 1.0;

	// printf("%e %e %e %e\n", p0, q0, p1, p2);
	// printf("%e %e %e\n", q1, q2, p3);

	double r0, r1, r2;
	// (p0 + p1 + p2 + q0) + (q1 + q2 + p3)
	p0 = qd::two_sum(p0, q0, r0);
	p1 = qd::two_sum(p1, p2, r1);
	p0 = qd::two_sum(p0, p1, r2);

	q0 = r0 + r1 + r2 + q1 + q2 + p3;
	
	p0 = qd::quick_two_sum(p0, q0, q0);

	return dd_real(p0, q0);
}

dd_real mult_sub1_qd(const dd_real &a, const dd_real &b){
	return to_dd_real(
			qd_real(a) * b - 1.0);
}

qd_real rsqrt_qd1(const dd_real &x){
	double y0 = 1.0 / sqrt(to_double(x));
	dd_real h = mult_sub1(x,  dd_real::sqr(y0));
	// dd_real c = dd_real(-0.5, (3./8.) * to_double(h));
	dd_real c = -0.5 + (3./8.)*h;

	qd_real y1 = y0 + qd_real(y0 * (h * c));

	return y1;
}
qd_real rsqrt_qd2(const dd_real &x){
	return 1.0 / sqrt(qd_real(x));
}

dd_real rsqrt_dd1(const dd_real &x){
	double y0 = 1.0 / sqrt(to_double(x));
	dd_real h = mult_sub1(x,  dd_real::sqr(y0));
	// dd_real c = dd_real(-0.5, (3./8.) * to_double(h));
	dd_real c = -0.5 + (3./8.)*h;

	dd_real y1 = y0 + (y0 * (h * c));

	return y1;
}
dd_real rsqrt_dd2(const dd_real &x){
	return to_dd_real(
			1.0 / sqrt(qd_real(x)));
}

int main(){
	using std::cout;
	using std::endl;

	srand(20170728);

	for(int i=0; i<10; i++){

		dd_real x  = ddrand();
#if 0
		double  y  = 1.0 / std::sqrt(to_double(x));
		dd_real y2 = dd_real::sqr(y);
		dd_real h1 = mult_sub1   (x, y2);
		dd_real h2 = mult_sub1_qd(x, y2);

		cout.precision(32);
		cout << h1 << endl;
		cout << h2 << endl;
		printf("%A %A\n", h1.x[0], h1.x[1]);
		printf("%A %A\n", h2.x[0], h2.x[1]);
		puts("");
#else
		dd_real y1 = rsqrt_dd1(x);
		dd_real y2 = rsqrt_dd2(x);
		cout.precision(32);
		cout << y1 << endl;
		cout << y2 << endl;
		printf("%A %A\n", y1.x[0], y1.x[1]);
		printf("%A %A\n", y2.x[0], y2.x[1]);
		puts("");
#endif
	}

	for(int i=0; i<1000000; i++){
		dd_real x  = ddrand();
		dd_real y1 = rsqrt_dd1(x);
		dd_real y2 = rsqrt_dd2(x);

		assert(y1 == y2);
	}

	return 0;
}
