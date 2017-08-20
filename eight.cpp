// eight-shaped 3-body problerm
#include <cstdio>
#include "vector3.h"

#include <qd/dd_real.h>
#include <qd/dd_inline.h>

inline dd_real rsqrt(const dd_real &x){
	double y_app = 1.0 / sqrt(to_double(x));
	dd_real x2 = mul_pwr2(x, 0.5);
	return y_app * (dd_real(1.5) - x2 * sqr(y_app));
}

typedef vector3<dd_real> qvec3;

template <> template <>
vector3<dd_real>::operator dvec3() const
{
	return dvec3(to_double(x), to_double(y), to_double(z));
}

#if defined HERMITE_FOURTH
#  include "hermite4.h"
#  include "test2.h"
#  define calc_force_on_i calc_force_on_i_p2
#elif defined HERMITE_SIXTH
#  include "hermite6.h"
#  include "test3.h"
#  define calc_force_on_i calc_force_on_i_p3
#elif defined HERMITE_EIGTH
#  include "hermite8.h"
#  include "test4.h"
#  define calc_force_on_i calc_force_on_i_p4
#elif defined HERMITE_TENTH
#  include "hermite10.h"
#  include "test5.h"
#  define calc_force_on_i calc_force_on_i_p5
#  define STABLE_INTERPOLATE
#elif defined HERMITE_TWELVETH
#  include "hermite12.h"
#  include "test6.h"
#  define calc_force_on_i calc_force_on_i_p6
#  define STABLE_INTERPOLATE
#elif defined HERMITE_FOURTEENTH
#  include "hermite14.h"
#  include "test7.h"
#  define calc_force_on_i calc_force_on_i_p7
#elif defined HERMITE_SIXTEENTH
#  include "hermite16.h"
#  include "test8.h"
#  define calc_force_on_i calc_force_on_i_p8
#else
#error
#endif
#include "nbodysystem.h"

int main(int ac, char **av){
	const dd_real dt = ac>1 ? 1./atoi(av[1]) : 1./1024.;
	const int npec   = ac>2 ? atoi(av[2])    : 1;
	const int norbit = ac>3 ? atoi(av[3])    : 10;

	NbodySystem sys;
	const int nbody = sys.nbody = 3;

	sys.ptcl .resize(nbody);
	sys.pred .resize(nbody);
	sys.force.resize(nbody);

	// position
	sys.ptcl[0].coord[0] = qvec3( 0.97000436, -0.24308753, 0.0);
	sys.ptcl[1].coord[0] = qvec3(-0.97000436, +0.24308753, 0.0);
	sys.ptcl[2].coord[0] = qvec3(0.0,          0.0,        0.0);
	// velocity
	sys.ptcl[0].coord[1] = qvec3(-0.93240737/(-2.0),-0.86473146/(-2.0),0.0);
	sys.ptcl[1].coord[1] = qvec3(-0.93240737/(-2.0),-0.86473146/(-2.0),0.0);
	sys.ptcl[2].coord[1] = qvec3(-0.93240737,       -0.86473146,       0.0);
	// mass
	sys.ptcl[0].mass = 1.0;
	sys.ptcl[1].mass = 1.0;
	sys.ptcl[2].mass = 1.0;

	sys.tsys = 0.0;

	sys.set_fixed_dt(dt);
	sys.init_force();
	const dd_real e0 = sys.calc_energy_from_ptcl();
	fprintf(stderr, "e0 : %24.16f\n", to_double(e0));

	const dd_real tperiod = 8. * atan(1.0);
	// const dd_real tperiod = dd_real::_2pi;

	dd_real err_max = 0.0;
	for(int i=0; ; i++){
		// for the first step, iterate many times
		const int n = i ? npec : npec+4;
		sys.tsys += dt;
		sys.predict_all();

		// P(EC)^n iteration
		for(int nn=0; nn<n; nn++){
			sys.calc_force_on_first_nact(nbody);
			sys.correct_and_feedback();
		}
		sys.calc_force_on_first_nact(nbody);
		sys.correct_and_commit();

		const dd_real e1 = sys.calc_energy_from_ptcl();
		const dd_real de = (e1-e0)/e0;
		err_max = std::max(err_max, fabs(de));
#if 1
		dvec3 r0 = sys.ptcl[0].coord[0];
		dvec3 r1 = sys.ptcl[1].coord[0];
		dvec3 r2 = sys.ptcl[2].coord[0];
		printf("%f %e  %f %f  %f %f  %f %f\n",
				to_double(sys.tsys), to_double(de),
				r0.x, r0.y, r1.x, r1.y, r2.x, r2.y);
		// exit(0);
#endif
		if(sys.tsys > norbit*tperiod) break;
	}
	const dd_real e1 = sys.calc_energy_from_ptcl();
	fprintf(stderr, "e1 : %24.16f\n", to_double(e1));
	printf("%e %e %d %d #grep\n", to_double(dt), to_double(err_max), norbit, npec);


	return 0;
}

