#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

struct NbodySystem{
	long   nbody;
	long   num_step,     num_bstep;
	long   num_step_tot, num_bstep_tot;

	typedef dd_real real_type;
	typedef qvec3   vect_type;

	real_type eps2;
	real_type tsys;
	real_type init_energy, prev_energy;
	real_type dtmax;
	real_type eta, eta_s;

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

		const real_type eps = 4.0 / nbody;
		this->eps2 = eps * eps;

		num_step = num_bstep = 0;
		num_step_tot = num_bstep_tot = 0;
	}
#endif

	real_type get_mass(const int i) const {
		return pred[i].mass;
	}

	template <int p> 
	vect_type get_coord(const int i) const {
		return pred[i].coord[p];
	}

	template <int p> 
	void set_force(const int i, const vect_type &f){
		force[i].force[p] = f;
	}

	template <typename PTCL> // Particle or Predictor
	real_type calc_energy(const std::vector<PTCL> &p){
		real_type pe  = 0.0;
		real_type ke2 = 0.0;

		for(int i=0; i<nbody; i++){
			real_type phi = 0.0;
			for(int j=i+1; j<nbody; j++){
				vect_type dr = p[j].coord[0] - p[i].coord[0];
				// real_type rinv = 1.0 / sqrt(eps2 + dr*dr);
				real_type rinv = rsqrt(eps2 + dr*dr);
				phi -= p[j].mass * rinv;
			}
			pe += p[i].mass * phi;
		}
		for(int i=0; i<nbody; i++){
			ke2 += p[i].mass * p[i].coord[1].norm2();
		}

		return pe + 0.5 * ke2;
	}

	real_type calc_energy_from_ptcl(){
		return calc_energy(ptcl);
	}
	real_type calc_energy_from_pred(){
		return calc_energy(pred);
	}

	void calc_force_on_first_nact(int nact){
#pragma omp parallel for
		for(int i=0; i<nact; i++){
			calc_force_on_i <NbodySystem, real_type, vect_type> 
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
	void set_fixed_dt(const real_type dt){
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
	void round_predictor(){
		for(int i=0; i<nbody; i++){
			for(int j=0; j<Particle::NFORCE; j++){
				pred[i].coord[j].x.x[1] = 0.0;
				pred[i].coord[j].y.x[1] = 0.0;
				pred[i].coord[j].z.x[1] = 0.0;
			}
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
	void integrate_pec_nth(const int n, const real_type dt){
		tsys += dt;
		predict_all();
		for(int iter=0; iter<n-1; iter++){
			calc_force_on_first_nact(nbody);
			correct_and_feedback();
		}
		calc_force_on_first_nact(nbody);
		correct_and_commit();
	}
};
