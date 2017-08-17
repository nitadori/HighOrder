#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

struct CmpPtcl_dt{
	bool operator()(const Particle &p1, const Particle &p2) const {
		return (p1.dt < p2.dt);
	}
};

struct NbodySystem{
	long   nbody;
	long   num_step,     num_bstep;
	long   num_step_tot, num_bstep_tot;

	typedef dd_real real_t;
	typedef qvec3   vect_t;

	real_t eps2;
	real_t tsys;
	real_t init_energy, prev_energy;
	real_t dtmax;
	real_t eta, eta_s;

	std::vector<Particle>  ptcl;
	std::vector<Predictor> pred;
	std::vector<Force>     force;

	void reset_counters(){
		init_energy = prev_energy;
		tsys = 0.0;
		num_step_tot = num_bstep_tot = 0;
		for(int i=0; i<nbody; i++) ptcl[i].tlast = tsys;
		puts("reset done!");
	}

#if 0
	void read_snapshot_masaki(const char *filename){
		int nread;
		FILE *fp = fopen(filename, "r");
		assert(fp);

		int snapid;
		long nbody;
		nread = fscanf(fp, "%d %ld %lf", &snapid, &nbody, &tsys);
		this->nbody = nbody;
		assert(3 == nread);
		fprintf(stderr, "read snapshot (form M.I.), n = %ld, t = %f\n", nbody, tsys);

		ptcl .resize(nbody);
		pred .resize(nbody);
		force.resize(nbody);

		for(int i=0; i<nbody; i++){
			long id;
			double mass, pot;
			dvec3 pos, vel;
			nread = fscanf(fp, "%ld %lA %lA %lA %lA %lA %lA %lA %lA", 
					&id, &mass,
					&pos.x, &pos.y, &pos.z,
					&vel.x, &vel.y, &vel.z,
					&pot);
			assert(9 == nread);

			ptcl[i].id       = id;
			ptcl[i].mass     = mass;
			ptcl[i].coord[0] = qvec3(pos);
			ptcl[i].coord[1] = qvec3(vel);
		}
		fclose(fp);

		const real_t eps = 4.0 / nbody;
		this->eps2 = eps * eps;

		num_step = num_bstep = 0;
		num_step_tot = num_bstep_tot = 0;
	}
#endif

	real_t get_mass(const int i) const {
		return pred[i].mass;
	}

	template <int p> 
	vect_t get_coord(const int i) const {
		return pred[i].coord[p];
	}

	template <int p> 
	void set_force(const int i, const vect_t &f){
		force[i].force[p] = f;
	}

	template <typename PTCL> // Particle or Predictor
	real_t calc_energy(const std::vector<PTCL> &p){
		real_t pe  = 0.0;
		real_t ke2 = 0.0;

		for(int i=0; i<nbody; i++){
			real_t phi = 0.0;
			for(int j=i+1; j<nbody; j++){
				vect_t dr = p[j].coord[0] - p[i].coord[0];
				// real_t rinv = 1.0 / sqrt(eps2 + dr*dr);
				real_t rinv = rsqrt(eps2 + dr*dr);
				phi -= p[j].mass * rinv;
			}
			pe += p[i].mass * phi;
		}
		for(int i=0; i<nbody; i++){
			ke2 += p[i].mass * p[i].coord[1].norm2();
		}

		return pe + 0.5 * ke2;
	}

	real_t calc_energy_from_ptcl(){
		return calc_energy(ptcl);
	}
	real_t calc_energy_from_pred(){
		return calc_energy(pred);
	}

	void calc_force_on_first_nact(int nact){
#pragma omp parallel for
		for(int i=0; i<nact; i++){
			calc_force_on_i <NbodySystem, real_t, vect_t> 
				(*this, nbody, i, this->eps2);
		}
	}
	void init_force(){
		const int niter = (Particle::NFORCE+1)/2;
		for(int n=0; n<niter; n++){
			for(int i=0; i<nbody; i++){
				pred[i].zero_predict(ptcl[i]);
			}
			calc_force_on_first_nact(nbody);
			for(int i=0; i<nbody; i++){
				force[i].init_assign(ptcl[i]);
			}
		}
	}
	void set_fixed_dt(const real_t dt){
		for(int i=0; i<nbody; i++){
			ptcl[i].tlast = tsys;
			ptcl[i].dt = dt;
		}
	}
	void predict_all(){
		for(int i=0; i<nbody; i++){
			pred[i].predict(tsys, ptcl[i]);
		}
	}
	void correct_and_feedback(){
		for(int i=0; i<nbody; i++){
			Corrector corr;
			corr.correct    (ptcl[i], force[i]);
			corr.feedback   (pred[i]);
		}
	}
	void correct_and_commit(){
		for(int i=0; i<nbody; i++){
			Corrector corr;
			corr.correct    (ptcl[i], force[i]);
			corr.interpolate(ptcl[i], force[i]);
#if 0
			corr.taylor_test(ptcl[i]);
			puts("");
			// exit(0);
#endif
			corr.commit     (ptcl[i]);
		}
	}
	void integrate_pec_nth(const int n, const real_t dt){
		tsys += dt;
		predict_all();
		for(int iter=0; iter<n-1; iter++){
			calc_force_on_first_nact(nbody);
			correct_and_feedback();
		}
		calc_force_on_first_nact(nbody);
		correct_and_commit();
	}
	// METHODS FOR INDIVIDUAL TIMESTEP SCHEME
#if 0
	int count_nact(const real_t tnext) const {
		int nact;
		for(nact=0; nact<nbody; nact++){
			if(tnext != ptcl[nact].tlast + ptcl[nact].dt)
				break;
		}
		return nact;
	}
	real_t calc_dtlim(const real_t tnext) const {
		real_t dtlim = dtmax;
		real_t s = tnext / dtmax;
		while(s != real_t(int(s))){
			s *= 2.0;
			dtlim *= 0.5;
			assert(dtlim >= 1.0/(1LL<<32));
		}
		return dtlim;
	}
	void sort_ptcl(const int nact){
		std::sort(&ptcl[0], &ptcl[nact], CmpPtcl_dt());
	}


	__attribute__((noinline))
	void integrate_one_block(){
		const real_t tnext = ptcl[0].tlast + ptcl[0].dt;
		const real_t dtlim = calc_dtlim(tnext);
		const int    nact  = count_nact(tnext);
#if 0
		printf("t = %f, nact = %6d, dtlim = %A\n", tsys, nact, dtlim);
#endif
		tsys = tnext;
		predict_all();

		calc_force_on_first_nact(nact);
#ifdef LOCAL_PECEC
		for(int i=0; i<nact; i++){
			Corrector corr;
			corr.correct    (ptcl[i], force[i]);
			corr.feedback   (pred[i]);
		}
		calc_force_on_first_nact(nact);
#endif

		for(int i=0; i<nact; i++){
			Corrector corr;
#ifdef STABLE_INTERPOLATE
			Corrector corr_stab;
			corr_stab.interpolate_stab(ptcl[i], force[i]);
#endif
			corr.correct    (ptcl[i], force[i]);
			corr.interpolate(ptcl[i], force[i]);
			corr.commit     (ptcl[i]);

			ptcl[i].calc_dt(eta, eta_s, dtlim);
#ifdef STABLE_INTERPOLATE
			corr_stab.overwrite_stab(ptcl[i]);
#endif
		}

		sort_ptcl(nact);
		// sort_ptcl(nact, dtlim);

		num_step += nact;
		num_bstep++;
	}
	void integrate_one_dtmax(){
		const real_t tt = tsys + dtmax;
		while(tsys < tt){
			integrate_one_block();
		}
		assert(tsys == tt);


		const real_t energy = calc_energy(ptcl);
		const real_t de_glo =((init_energy - energy) / init_energy);
		const real_t de_loc =((prev_energy - energy) / init_energy);
		assert(de_glo < 1.0);
		num_step_tot  += num_step;
		num_bstep_tot += num_bstep;
		const double nact_loc = double(num_step) / double(num_bstep);
		const double nact_glo = double(num_step_tot) / double(num_bstep_tot);

		fprintf(stderr, "t = %f\n", tsys);
		fprintf(stderr, " steps: %ld %ld %ld %ld\n", num_bstep, num_step, num_bstep_tot, num_step_tot);
		fprintf(stderr, " nact : %f %f\n", nact_loc, nact_glo);
		fprintf(stderr, " de(local/global) : %+e %+e\n", to_double(de_loc), to_double(de_glo));

		FILE *fp = stdout;
		fprintf(fp, "%f %6ld %6ld %10ld %10ld  %8.2f %8.2f   %+e  %+e\n", 
				tsys, 
				num_bstep, num_step, 
				num_bstep_tot, num_step_tot, 
				nact_loc, nact_glo,
				de_loc, de_glo); 
		fflush(fp);

		prev_energy = energy;
		num_step  = 0;
		num_bstep = 0;
	}
#endif
};
