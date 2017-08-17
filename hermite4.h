#include "vector3.h"
#include "taylor.h"

struct Particle{
	enum
	{
		NFORCE = 2,
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
	qvec3 coord[2];

	void predict(const dd_real tnext, const Particle &p){
		const dd_real dt = tnext - p.tlast;
		mass = p.mass;
		coord[0] = taylor<3>(dt, &p.coord[0]);
		coord[1] = taylor<2>(dt, &p.coord[1]);
	}
	void zero_predict(const Particle &p){
		mass     = p.mass;
		coord[0] = p.coord[0];
		coord[1] = p.coord[1];
	}
};

struct Force{
	qvec3 force[2]; // acc, jerk
	void init_assign(Particle &p){
		for(int m=0; m<Particle::NFORCE; m++){
			p.coord[2+m] = force[m];
		}
	}
};

struct Corrector{
	qvec3 pos;
	qvec3 vel;
	qvec3 force[2];
	qvec3 fintp[2];

	void correct(const Particle &p, const Force &f){
		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;

		const dd_real d0 = c0;
		const dd_real d1 = dd_real::div(1., 3.)  * c1;

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// correct velocity
		const qvec3 f0pl = d0 * (fr[0] + fl[0]);
		const qvec3 f1mn = d1 * (fr[1] - fl[1]);

		this->vel = p.coord[1] + p.dt * (f0pl - (f1mn));

		// correct position
		const qvec3 v0pl = d0 * (this->vel + p.coord[1]);
		const qvec3 v1mn = d1 * (fr[0] - fl[0]);

		this->pos = p.coord[0] + p.dt * (v0pl - (v1mn));

		// copy forces
		for(int i=0; i<2; i++) force[i] = f.force[i];
	}	
	void interpolate(const Particle &p, const Force &f){
		qvec3 fpl[2];
		qvec3 fmn[2];
		qvec3 fmid[4];
		// qvec3 fleft [4];
		qvec3 fright[4];

		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// f+, f-
		{
			fpl[0] = c0 * (fr[0] + fl[0]);
			fmn[0] = c0 * (fr[0] - fl[0]);
			fpl[1] = c1 * (fr[1] + fl[1]);
			fmn[1] = c1 * (fr[1] - fl[1]);
		}
		// fmid
		{
			fmid[2] = 1./2. * fmn[1];
			fmid[3] = 1./2. * (-fmn[0] + fpl[1]);
		}
		fright[2] = fmid[2] +  3.*fmid[3];
		fright[3] =               fmid[3];

		// rescale
		{
			const dd_real hi = 1.0 / h;
			const dd_real s1 = hi;
			const dd_real s2 = s1 * (2.*hi);
			const dd_real s3 = s2 * (3.*hi);

			fintp[0] = s2 * fright[2];
			fintp[1] = s3 * fright[3];
		}
	}	
	void commit(Particle &p){
		p.coord[0] = pos;
		p.coord[1] = vel;

		p.coord[2] = force[0];
		p.coord[3] = force[1];

		p.coord[4] = fintp[0];
		p.coord[5] = fintp[1];

		p.tlast += p.dt;
	}
	void feedback(Predictor &pr){
		pr.coord[0] = pos;
		pr.coord[1] = vel;
	}
	void taylor_test(const Particle &p){
		dvec3 df0 = p.coord[2] - taylor<3>(-p.dt, &force[0]);
		dvec3 df1 = p.coord[3] - taylor<2>(-p.dt, &force[1]);
		printf("[0]: %e %e %e\n", df0.x, df0.y, df0.z);
		printf("[1]: %e %e %e\n", df1.x, df1.y, df1.z);
	}
};
