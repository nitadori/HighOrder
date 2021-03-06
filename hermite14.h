#include "vector3.h"
#include "taylor.h"

struct Particle{
	enum
	{
		NFORCE = 7,
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
	qvec3 coord[7];

	void predict(const dd_real tnext, const Particle &p){
		const dd_real dt = tnext - p.tlast;
		mass = p.mass;
		coord[0] = taylor<13>(dt, &p.coord[0]);
		coord[1] = taylor<12>(dt, &p.coord[1]);
		coord[2] = taylor<11>(dt, &p.coord[2]);
		coord[3] = taylor<10>(dt, &p.coord[3]);
		coord[4] = taylor< 9>(dt, &p.coord[4]);
		coord[5] = taylor< 8>(dt, &p.coord[5]);
		coord[6] = taylor< 7>(dt, &p.coord[6]);
	}
	void zero_predict(const Particle &p){
		mass     = p.mass;
		coord[0] = p.coord[0];
		coord[1] = p.coord[1];
		coord[2] = p.coord[2];
		coord[3] = p.coord[3];
		coord[4] = p.coord[4];
		coord[5] = p.coord[5];
		coord[6] = p.coord[6];
	}
};

struct Force{
	qvec3 force[7]; // acc, jerk, snap, crackle, pop, pdot, pdot2
	void init_assign(Particle &p){
		for(int m=0; m<Particle::NFORCE; m++){
			p.coord[2+m] = force[m];
		}
	}
};

struct Corrector{
	qvec3 pos;
	qvec3 vel;
	qvec3 force[7];
	qvec3 fintp[7];

	void correct(const Particle &p, const Force &f){
		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * (dd_real::div(1., 2.) * h);
		const dd_real c3 = c2 * (dd_real::div(1., 3.) * h);
		const dd_real c4 = c3 * (dd_real::div(1., 4.) * h);
		const dd_real c5 = c4 * (dd_real::div(1., 5.) * h);
		const dd_real c6 = c5 * (dd_real::div(1., 6.) * h);

#if 1
		const dd_real d0 = c0;
		const dd_real d1 = dd_real::div( 6.,   13.) * c1;
		const dd_real d2 = dd_real::div(10.,   39.) * c2;
		const dd_real d3 = dd_real::div(20.,  143.) * c3;
		const dd_real d4 = dd_real::div(48.,  715.) * c4;
		const dd_real d5 = dd_real::div(32., 1287.) * c5;
		const dd_real d6 = dd_real::div(16., 3003.) * c6;
#else
#endif

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// correct velocity
		const qvec3 f0pl = d0 * (fr[0] + fl[0]);
		const qvec3 f1mn = d1 * (fr[1] - fl[1]);
		const qvec3 f2pl = d2 * (fr[2] + fl[2]);
		const qvec3 f3mn = d3 * (fr[3] - fl[3]);
		const qvec3 f4pl = d4 * (fr[4] + fl[4]);
		const qvec3 f5mn = d5 * (fr[5] - fl[5]);
		const qvec3 f6pl = d6 * (fr[6] + fl[6]);

		this->vel = p.coord[1] + p.dt * (f0pl - (f1mn - (f2pl - (f3mn - (f4pl - (f5mn - f6pl))))));

		// correct position
		const qvec3 v0pl = d0 * (this->vel + p.coord[1]);
		const qvec3 v1mn = d1 * (fr[0] - fl[0]);
		const qvec3 v2pl = d2 * (fr[1] + fl[1]);
		const qvec3 v3mn = d3 * (fr[2] - fl[2]);
		const qvec3 v4pl = d4 * (fr[3] + fl[3]);
		const qvec3 v5mn = d5 * (fr[4] - fl[4]);
		const qvec3 v6pl = d6 * (fr[5] + fl[5]);

		this->pos = p.coord[0] + p.dt * (v0pl - (v1mn - (v2pl - (v3mn - (v4pl - (v5mn - v6pl))))));

		// copy forces
		for(int i=0; i<7; i++) force[i] = f.force[i];
	}

	void interpolate(const Particle &p, const Force &f){
		qvec3 fpl[7];
		qvec3 fmn[7];
		qvec3 fmid[14];
		// qvec3 fleft [14];
		qvec3 fright[14];

		const dd_real h = 0.5 * p.dt;
		const dd_real c0 = 0.5;
		const dd_real c1 = c0 * h;
		const dd_real c2 = c1 * (dd_real::div(1., 2.) * h);
		const dd_real c3 = c2 * (dd_real::div(1., 3.) * h);
		const dd_real c4 = c3 * (dd_real::div(1., 4.) * h);
		const dd_real c5 = c4 * (dd_real::div(1., 5.) * h);
		const dd_real c6 = c5 * (dd_real::div(1., 6.) * h);

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
			fpl[5] = c5 * (fr[5] + fl[5]);
			fmn[5] = c5 * (fr[5] - fl[5]);
			fpl[6] = c6 * (fr[6] + fl[6]);
			fmn[6] = c6 * (fr[6] - fl[6]);
		}
#if 0
		// fmid
		{
			qvec3 evn[7], odd[7];

			evn[1] =          fmn[1] - 2.*fpl[2] + 3.*fmn[3] - 4.*fpl[4] +  5.*fmn[5] -  6.*fpl[6];
			evn[3] =                                  fmn[3] - 4.*fpl[4] + 10.*fmn[5] - 20.*fpl[6];
			evn[5] =                                                           fmn[5] -  6.*fpl[6];

			odd[0] = fmn[0] - fpl[1] +    fmn[2] -    fpl[3] +    fmn[4] -     fpl[5] +     fmn[6];
			odd[2] =                      fmn[2] - 3.*fpl[3] + 6.*fmn[4] - 10.*fpl[5] + 15.*fmn[6];
			odd[4] =                                              fmn[4] -  5.*fpl[5] + 15.*fmn[6];
			odd[6] =                                                                        fmn[6];

			fmid[ 8] = (1./1024.) * (-495.*evn[1] + 189.*evn[3] - 175.*evn[5]);
			fmid[10] = (1./1024.) * ( 154.*evn[1] -  54.*evn[3] +  42.*evn[5]);
			fmid[12] = (1./1024.) * ( -21.*evn[1] +   7.*evn[3] -   5.*evn[5]);

			fmid[ 7] = (1./1024.) * (-8580.*odd[0] +  924.*odd[2] - 420.*odd[4] + 700.*odd[6]);
			fmid[ 9] = (1./1024.) * ( 5005.*odd[0] -  495.*odd[2] + 189.*odd[4] - 175.*odd[6]);
			fmid[11] = (1./1024.) * (-1638.*odd[0] +  154.*odd[2] -  54.*odd[4] +  42.*odd[6]);
			fmid[13] = (1./1024.) * (  231.*odd[0] -   21.*odd[2] +   7.*odd[4] -   5.*odd[6]);
		}
#else
		// even
		{
			qvec3 tmp  = fpl[2] - (1./2.) * fmn[1];

			// lower triangle
			qvec3 evn0 = (1./16.) * ((5./4.)  * tmp + (-3./2.) * fmn[3] + fpl[4]);
			qvec3 evn1 = (1./32.) * ((-7./4.) * tmp + (9./4)   * fmn[3] + (-2.0)  * fpl[4] + fmn[5]);
			qvec3 evn2 = (1./64.) * ((21./8.) * tmp + (-7./2.) * fmn[3] + (7./2.) * fpl[4] + (-5./2.) * fmn[5] + fpl[6]);

			// upper triangle
			fmid[8]  = evn0 - 5.0 * evn1 + 15.0 * evn2;
			fmid[10] =              evn1 -  6.0 * evn2;
			fmid[12] =                            evn2;
		}
		// odd
		{
			qvec3 tmp = fpl[1] - fmn[0];

			// lower triangle
			qvec3 odd0 = (1./8.)  * ((5./2.)     * tmp + (-2.0)    * fmn[2] +             fpl[3]);
			qvec3 odd1 = (1./16.) * ((-35./8.)   * tmp + (15./4.)  * fmn[2] + (-5./2.)  * fpl[3] +          fmn[4]);
			qvec3 odd2 = (1./32.) * ((63./8)     * tmp + (-7.0)    * fmn[2] + (21./4)   * fpl[3] + (-3.0) * fmn[4] + fpl[5]);
			qvec3 odd3 = (1./64.) * ((-231./16.) * tmp + (105./8.) * fmn[2] + (-21./2.) * fpl[3] + (7.0) * fmn[4] + (-7./2.) * fpl[5] + fmn[6]);

			// upper triangle
			fmid[7]  = odd0 - 4.0 * odd1 + 10.0 * odd2 - 20.0 * odd3;
			fmid[9]  =              odd1 - 5.0  * odd2 + 15.0 * odd3;
			fmid[11] =                            odd2 -  6.0 * odd3;
			fmid[13] =                                          odd3;
		}
#endif
		fright[ 7] = fmid[7] +  8.*fmid[8] +  36.*fmid[9] + 120.*fmid[10] + 330.*fmid[11] + 792.*fmid[12] + 1716.*fmid[13];
		fright[ 8] =               fmid[8] +   9.*fmid[9] +  45.*fmid[10] + 165.*fmid[11] + 495.*fmid[12] + 1287.*fmid[13];
		fright[ 9] =                              fmid[9] +  10.*fmid[10] +  55.*fmid[11] + 220.*fmid[12] +  715.*fmid[13];
		fright[10] =                                             fmid[10] +  11.*fmid[11] +  66.*fmid[12] +  286.*fmid[13];
		fright[11] =                                                             fmid[11] +  12.*fmid[12] +   78.*fmid[13];
		fright[12] =                                                                             fmid[12] +   13.*fmid[13];
		fright[13] =                                                                                              fmid[13];
#if 0 // 12th-order interpolation for stability
#endif
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
			const dd_real s10= s9 *(10.*hi);
			const dd_real s11= s10*(11.*hi);
			const dd_real s12= s11*(12.*hi);
			const dd_real s13= s12*(13.*hi);

			fintp[0] = s7 * fright[7];
			fintp[1] = s8 * fright[8];
			fintp[2] = s9 * fright[9];
			fintp[3] = s10* fright[10];
			fintp[4] = s11* fright[11];
			fintp[5] = s12* fright[12];
			fintp[6] = s13* fright[13];
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
		p.coord[7] = force[5];
		p.coord[8] = force[6];

		p.coord[ 9] = fintp[0];
		p.coord[10] = fintp[1];
		p.coord[11] = fintp[2];
		p.coord[12] = fintp[3];
		p.coord[13] = fintp[4];
		p.coord[14] = fintp[5];
		p.coord[15] = fintp[6];

		p.tlast += p.dt;
	}
	void feedback(Predictor &pr){
		pr.coord[0] = pos;
		pr.coord[1] = vel;
		pr.coord[2] = force[0]; // acc
		pr.coord[3] = force[1]; // jrk
		pr.coord[4] = force[2]; // snp
		pr.coord[5] = force[3]; // crk
		pr.coord[6] = force[4]; // pop
	}
	void taylor_test(const Particle &p){
		dvec3 df0 = p.coord[2] - taylor<13>(-p.dt, &force[0]);
		dvec3 df1 = p.coord[3] - taylor<12>(-p.dt, &force[1]);
		dvec3 df2 = p.coord[4] - taylor<11>(-p.dt, &force[2]);
		dvec3 df3 = p.coord[5] - taylor<10>(-p.dt, &force[3]);
		dvec3 df4 = p.coord[6] - taylor< 9>(-p.dt, &force[4]);
		dvec3 df5 = p.coord[7] - taylor< 8>(-p.dt, &force[5]);
		dvec3 df6 = p.coord[8] - taylor< 7>(-p.dt, &force[6]);
		printf("[0]: %e %e %e\n", df0.x, df0.y, df0.z);
		printf("[1]: %e %e %e\n", df1.x, df1.y, df1.z);
		printf("[2]: %e %e %e\n", df2.x, df2.y, df2.z);
		printf("[3]: %e %e %e\n", df3.x, df3.y, df3.z);
		printf("[4]: %e %e %e\n", df4.x, df4.y, df4.z);
		printf("[5]: %e %e %e\n", df5.x, df5.y, df5.z);
		printf("[6]: %e %e %e\n", df6.x, df6.y, df6.z);
	}
#if 0
	void interpolate_stab(const Particle &p, const Force &f){
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
	// 10th-order interpolation for stability
		{
			qvec3 evn[5], odd[5];
			evn[1] = fmn[1] - 2.*fpl[2] + 3.*fmn[3] - 4.*fpl[4];
			evn[3] =                         fmn[3] - 4.*fpl[4];

			odd[0] = fmn[0] - fpl[1] +  fmn[2] -    fpl[3] +    fmn[4];
			odd[2] =                    fmn[2] - 3.*fpl[3] + 6.*fmn[4];
			odd[4] =                                            fmn[4];

			fmid[6] = (1./32. ) * (  7.*evn[1] -  5.*evn[3]);
			fmid[8] = (1./128.) * ( -5.*evn[1] +  3.*evn[3]);

			fmid[5] = (1./64. ) * ( 189.*odd[0] - 35.*odd[2] + 45.*odd[4]);
			fmid[7] = (1./32. ) * ( -45.*odd[0] +  7.*odd[2] -  5.*odd[4]);
			fmid[9] = (1./128.) * (  35.*odd[0] -  5.*odd[2] +  3.*odd[4]);
		}
		fright[6] = fmid[6] +  7.*fmid[7] + 28.*fmid[8] +  84.*fmid[9];
		fright[7] =               fmid[7] +  8.*fmid[8] +  36.*fmid[9];
		fright[8] =                             fmid[8] +   9.*fmid[9];
		fright[9] =                                            fmid[9];
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

			fintp[0] = s6 * fright[6];
			fintp[1] = s7 * fright[7];
			fintp[2] = s8 * fright[8];
			fintp[3] = s9 * fright[9];
		}
	}
	void overwrite_stab(Particle &p){
		p.coord[ 8] = fintp[0];
		p.coord[ 9] = fintp[1];
		p.coord[10] = fintp[2];
		p.coord[11] = fintp[3];
	}
#endif
};

