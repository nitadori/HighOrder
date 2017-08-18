#include "vector3.h"
#include "taylor.h"

struct Particle{
	enum
	{
		NFORCE = 3,
		ORDER  = 3*NFORCE,
	};
	unsigned long id;
	dd_real dt;
	dd_real tlast;
	dd_real mass;

	qvec3 coord[2 + ORDER]; // pos-vel, forces, interpolates
};

struct Predictor{
	dd_real mass;
	qvec3 coord[3];

	void predict(const dd_real tnext, const Particle &p){
		const dd_real dt = tnext - p.tlast;
		mass = p.mass;
		coord[0] = taylor<5>(dt, &p.coord[0]);
		coord[1] = taylor<4>(dt, &p.coord[1]);
		coord[2] = taylor<3>(dt, &p.coord[2]);
	}
	void zero_predict(const Particle &p){
		mass     = p.mass;
		coord[0] = p.coord[0];
		coord[1] = p.coord[1];
		coord[2] = p.coord[2];
	}
};

struct Force{
	qvec3 force[3]; // acc, jerk, snap
	void init_assign(Particle &p){
		for(int m=0; m<Particle::NFORCE; m++){
			p.coord[2+m] = force[m];
		}
	}
};

struct Corrector{
	qvec3 pos;
	qvec3 vel;
	qvec3 force[3];
	qvec3 fintp[3];

	void correct(const Particle &p, const Force &f){
		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * ((1./2.) * h);

		const dd_real d0 = c0;
		const dd_real d1 = dd_real::div(2., 5.)  * c1;
		const dd_real d2 = dd_real::div(2., 15.) * c2;

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// correct velocity
		const qvec3 f0pl = d0 * (fr[0] + fl[0]);
		const qvec3 f1mn = d1 * (fr[1] - fl[1]);
		const qvec3 f2pl = d2 * (fr[2] + fl[2]);

		this->vel = p.coord[1] + p.dt * (f0pl - (f1mn - (f2pl)));

		// correct position
		const qvec3 v0pl = d0 * (this->vel + p.coord[1]);
		const qvec3 v1mn = d1 * (fr[0] - fl[0]);
		const qvec3 v2pl = d2 * (fr[1] + fl[1]);

		this->pos = p.coord[0] + p.dt * (v0pl - (v1mn - (v2pl)));

		// copy forces
		for(int i=0; i<3; i++) force[i] = f.force[i];
	}	
	void interpolate(const Particle &p, const Force &f){
		qvec3 fpl[3];
		qvec3 fmn[3];
		qvec3 fmid[6];
		// qvec3 fleft [6];
		qvec3 fright[6];

		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * ((1./2.) * h);

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
		}
		// fmid
		{
			fmid[3] = 1./8. * (-10. * fmn[0] + 10.* fpl[1] - 4. * fmn[2]);
			fmid[5] = 1./8. * (  3. * fmn[0] -  3.* fpl[1] + 2. * fmn[2]);
			fmid[4] = 1./8. *  (-fmn[1] + 2. * fpl[2]);
		}
		fright[3] = fmid[3] +  4.*fmid[4] + 10.*fmid[5];
		fright[4] =               fmid[4] +  5.*fmid[5];
		fright[5] =                             fmid[5];

		// rescale
		{
			const dd_real hi = 1.0 / h;
			const dd_real s1 = hi;
			const dd_real s2 = s1 * (2.*hi);
			const dd_real s3 = s2 * (3.*hi);
			const dd_real s4 = s3 * (4.*hi);
			const dd_real s5 = s4 * (5.*hi);

			fintp[0] = s3 * fright[3];
			fintp[1] = s4 * fright[4];
			fintp[2] = s5 * fright[5];
		}
	}	
	void commit(Particle &p){
		p.coord[0] = pos;
		p.coord[1] = vel;

		p.coord[2] = force[0];
		p.coord[3] = force[1];
		p.coord[4] = force[2];

		p.coord[5] = fintp[0];
		p.coord[6] = fintp[1];
		p.coord[7] = fintp[2];

		p.tlast += p.dt;
	}
	void feedback(Predictor &pr){
		pr.coord[0] = pos;
		pr.coord[1] = vel;
		pr.coord[2] = force[0]; // acc
	}
	void taylor_test(const Particle &p){
		dvec3 df0 = p.coord[2] - taylor<5>(-p.dt, &force[0]);
		dvec3 df1 = p.coord[3] - taylor<4>(-p.dt, &force[1]);
		dvec3 df2 = p.coord[4] - taylor<3>(-p.dt, &force[2]);
		printf("[0]: %e %e %e\n", df0.x, df0.y, df0.z);
		printf("[1]: %e %e %e\n", df1.x, df1.y, df1.z);
		printf("[2]: %e %e %e\n", df2.x, df2.y, df2.z);
	}
};
