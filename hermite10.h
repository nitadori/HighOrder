#include "vector3.h"
#include "taylor.h"

struct Particle{
	enum
	{
		NFORCE = 5,
		ORDER  = 2*NFORCE,
	};
	unsigned long id;
	dd_real dt;
	dd_real tlast;
	dd_real mass;

	qvec3 coord[2 + ORDER]; // pos-vel, forces, interpolates
};

struct Predictor{
	dd_real mass;
	qvec3 coord[5];

	void predict(const dd_real tnext, const Particle &p){
		const dd_real dt = tnext - p.tlast;
		mass = p.mass;

		coord[0] = taylor<9>(dt, &p.coord[0]);
		coord[1] = taylor<8>(dt, &p.coord[1]);
		coord[2] = taylor<7>(dt, &p.coord[2]);
		coord[3] = taylor<6>(dt, &p.coord[3]);
		coord[4] = taylor<5>(dt, &p.coord[4]);
	}
	void zero_predict(const Particle &p){
		mass     = p.mass;
		coord[0] = p.coord[0];
		coord[1] = p.coord[1];
		coord[2] = p.coord[2];
		coord[3] = p.coord[3];
		coord[4] = p.coord[4];
	}
};

struct Force{
	qvec3 force[5]; // acc, jerk, snap, crackle, pop
	void init_assign(Particle &p){
		for(int m=0; m<Particle::NFORCE; m++){
			p.coord[2+m] = force[m];
		}
	}
};

struct Corrector{
	qvec3 pos;
	qvec3 vel;
	qvec3 force[5];
	qvec3 fintp[5];

	void correct(const Particle &p, const Force &f){
		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * (dd_real::div(1., 2.) * h);
		const dd_real c3 = c2 * (dd_real::div(1., 3.) * h);
		const dd_real c4 = c3 * (dd_real::div(1., 4.) * h);

		const dd_real d0 = c0;
		const dd_real d1 = dd_real::div(4., 9.)   * c1;
		const dd_real d2 = dd_real::div(6., 27.)  * c2;
		const dd_real d3 = dd_real::div(2., 21.)  * c3;
		const dd_real d4 = dd_real::div(8., 315.) * c4;

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// correct velocity
		const qvec3 f0pl = d0 * (fr[0] + fl[0]);
		const qvec3 f1mn = d1 * (fr[1] - fl[1]);
		const qvec3 f2pl = d2 * (fr[2] + fl[2]);
		const qvec3 f3mn = d3 * (fr[3] - fl[3]);
		const qvec3 f4pl = d4 * (fr[4] + fl[4]);

		this->vel = p.coord[1] + p.dt * (f0pl - (f1mn - (f2pl - (f3mn - f4pl))));

		// correct position
		const qvec3 v0pl = d0 * (this->vel + p.coord[1]);
		const qvec3 v1mn = d1 * (fr[0] - fl[0]);
		const qvec3 v2pl = d2 * (fr[1] + fl[1]);
		const qvec3 v3mn = d3 * (fr[2] - fl[2]);
		const qvec3 v4pl = d4 * (fr[3] + fl[3]);

		this->pos = p.coord[0] + p.dt * (v0pl - (v1mn - (v2pl - (v3mn - v4pl))));

		// copy forces
		for(int i=0; i<5; i++) force[i] = f.force[i];
	}
	void interpolate(const Particle &p, const Force &f){
		qvec3 fpl[5];
		qvec3 fmn[5];
		qvec3 fmid[10];
		qvec3 fright[10];

		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * (dd_real::div(1., 2.) * h);
		const dd_real c3 = c2 * (dd_real::div(1., 3.) * h);
		const dd_real c4 = c3 * (dd_real::div(1., 4.) * h);

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// f+, f-
		{
			fpl[0] = c0 * (fr[0] + fl[0]);
			fmn[0] = c0 * (fr[0] - fl[0]);
			fpl[1] = c1 * (fr[1] + fl[1]);
			fmn[1] = c1 * (fr[1] - fl[1]);
			fpl[2] = c2 * (fr[2] + fl[2]);
			fmn[2] = c2 * (fr[2] - fl[2]);
			fpl[3] = c3 * (fr[3] + fl[3]);
			fmn[3] = c3 * (fr[3] - fl[3]);
			fpl[4] = c4 * (fr[4] + fl[4]);
			fmn[4] = c4 * (fr[4] - fl[4]);
		}
		// fmid
		{
			qvec3 evn[5], odd[5];
			evn[1] =          fmn[1] - 2.*fpl[2] + 3.*fmn[3] - 4.*fpl[4];
			evn[3] =                                  fmn[3] - 4.*fpl[4];

			odd[0] = fmn[0] - fpl[1] +    fmn[2] -    fpl[3] +    fmn[4];
			odd[2] =                      fmn[2] - 3.*fpl[3] + 6.*fmn[4];
			odd[4] =                                              fmn[4];

			fmid[6] = (1./32. ) * (  7.*evn[1] -  5.*evn[3]);
			fmid[8] = (1./128.) * ( -5.*evn[1] +  3.*evn[3]);

			fmid[5] = (1./64. ) * ( 189.*odd[0] - 35.*odd[2] + 45.*odd[4]);
			fmid[7] = (1./32. ) * ( -45.*odd[0] +  7.*odd[2] -  5.*odd[4]);
			fmid[9] = (1./128.) * (  35.*odd[0] -  5.*odd[2] +  3.*odd[4]);
		}

		fright[5] = fmid[5] +  6.*fmid[6] + 21.*fmid[7] + 56.*fmid[8] + 126.*fmid[9];
		fright[6] =               fmid[6] +  7.*fmid[7] + 28.*fmid[8] +  84.*fmid[9];
		fright[7] =                             fmid[7] +  8.*fmid[8] +  36.*fmid[9];
		fright[8] =                                           fmid[8] +   9.*fmid[9];
		fright[9] =                                                          fmid[9];

		// rescale
		{
			const dd_real hi = 1.0 / h;
			const dd_real s1 = hi;
			const dd_real s2 = s1 * (2.*hi);
			const dd_real s3 = s2 * (3.*hi);
			const dd_real s4 = s3 * (4.*hi);
			const dd_real s5 = s4 * (5.*hi);
			const dd_real s6 = s5 * (6.*hi);
			const dd_real s7 = s6 * (7.*hi);
			const dd_real s8 = s7 * (8.*hi);
			const dd_real s9 = s8 * (9.*hi);

			fintp[0] = s5 * fright[5];
			fintp[1] = s6 * fright[6];
			fintp[2] = s7 * fright[7];
			fintp[3] = s8 * fright[8];
			fintp[4] = s9 * fright[9];
		}
	}	
	void commit(Particle &p){
		p.coord[0] = pos;
		p.coord[1] = vel;

		p.coord[2] = force[0];
		p.coord[3] = force[1];
		p.coord[4] = force[2];
		p.coord[5] = force[3];
		p.coord[6] = force[4];

		p.coord[ 7] = fintp[0];
		p.coord[ 8] = fintp[1];
		p.coord[ 9] = fintp[2];
		p.coord[10] = fintp[3];
		p.coord[11] = fintp[4];

		p.tlast += p.dt;
	}
	void feedback(Predictor &pr){
		pr.coord[0] = pos;
		pr.coord[1] = vel;
		pr.coord[2] = force[0]; // acc
		pr.coord[3] = force[1]; // jrk
		pr.coord[4] = force[2]; // snp
	}
	void taylor_test(const Particle &p){
		dvec3 df0 = p.coord[2] - taylor<9>(-p.dt, &force[0]);
		dvec3 df1 = p.coord[3] - taylor<8>(-p.dt, &force[1]);
		dvec3 df2 = p.coord[4] - taylor<7>(-p.dt, &force[2]);
		dvec3 df3 = p.coord[5] - taylor<6>(-p.dt, &force[3]);
		dvec3 df4 = p.coord[6] - taylor<5>(-p.dt, &force[4]);
		printf("[0]: %e %e %e\n", df0.x, df0.y, df0.z);
		printf("[1]: %e %e %e\n", df1.x, df1.y, df1.z);
		printf("[2]: %e %e %e\n", df2.x, df2.y, df2.z);
		printf("[3]: %e %e %e\n", df3.x, df3.y, df3.z);
		printf("[4]: %e %e %e\n", df4.x, df4.y, df4.z);
	}
	void interpolate_stab(const Particle &p, const Force &f){
		qvec3 fpl[4];
		qvec3 fmn[4];
		qvec3 fmid[8];
		qvec3 fright[8];

		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * (dd_real::div(1., 2.) * h);
		const dd_real c3 = c2 * (dd_real::div(1., 3.) * h);

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// f+, f-
		{
			fpl[0] = c0 * (fr[0] + fl[0]);
			fmn[0] = c0 * (fr[0] - fl[0]);
			fpl[1] = c1 * (fr[1] + fl[1]);
			fmn[1] = c1 * (fr[1] - fl[1]);
			fpl[2] = c2 * (fr[2] + fl[2]);
			fmn[2] = c2 * (fr[2] - fl[2]);
			fpl[3] = c3 * (fr[3] + fl[3]);
			fmn[3] = c3 * (fr[3] - fl[3]);
		}
	// use 8th-order interpolation for stability
		{
			qvec3 evn[4], odd[4];
			evn[1] = fmn[1] - 2.*fpl[2] + 3.*fmn[3];
			evn[3] =                         fmn[3];

			odd[0] = fmn[0] - fpl[1] + fmn[2] -    fpl[3];
			odd[2] =                   fmn[2] - 3.*fpl[3];

			fmid[4] = (1./16.) * (-5.*evn[1] + 9.*evn[3]);
			fmid[6] = (1./16.) * (    evn[1] -    evn[3]);

			fmid[5] = (1./16.) * ( 21.*odd[0] -  5.*odd[2]);
			fmid[7] = (1./16.) * ( -5.*odd[0] +     odd[2]);
		}
		fright[4] = fmid[4] +
		         5.*fmid[5] + 15.*fmid[6] + 35.*fmid[7];
		fright[5] = fmid[5] +  6.*fmid[6] + 21.*fmid[7];
		fright[6] =               fmid[6] +  7.*fmid[7];
		fright[7] =                             fmid[7];
		// rescale
		{
			const dd_real hi = 1.0 / h;
			const dd_real s1 = hi;
			const dd_real s2 = s1 * (2.*hi);
			const dd_real s3 = s2 * (3.*hi);
			const dd_real s4 = s3 * (4.*hi);
			const dd_real s5 = s4 * (5.*hi);
			const dd_real s6 = s5 * (6.*hi);
			const dd_real s7 = s6 * (7.*hi);

			force[3] = s4 * fright[4];
			fintp[0] = s5 * fright[5];
			fintp[1] = s6 * fright[6];
			fintp[2] = s7 * fright[7];
		}
	}
	void overwrite_stab(Particle &p){
#if 0
		p.coord[6] = force[3];
#endif
		p.coord[7] = fintp[0];
		p.coord[8] = fintp[1];
		p.coord[9] = fintp[2];
	}
};
