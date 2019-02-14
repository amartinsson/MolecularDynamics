#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP
#include <iostream>
#include <vector>
#include <math.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <sys/time.h>
#include <omp.h>

#include "Generator.hpp"
#include "Grid.hpp"
#include "InfiniteSwitchSimulatedTempering.hpp"
#include "NptGrid.hpp"
#include "Particle.hpp"
#include "SimulatedTempering.hpp"
#include "System.hpp"

using namespace::std;

 /******************************************************************************
                            Generic Base Class
  *****************************************************************************/
class Integrator
{
public:
    Integrator() {};
    // destructor
    ~Integrator() {};
    // must be implemented
    virtual void integrate(Molecule* molecule_pt) = 0;
};

/******************************************************************************
                           Langevin Base Class
 *****************************************************************************/
class Langevin : public Integrator
{
public:
    Langevin(const double& beta, const double& gamma, const double& gamma_rot,
             const double& o_step_size, System* system_pt, const int& seed);
    // destructor
    ~Langevin();
    // integrator
    void integrate(Molecule* molecule_pt);
    // set to integrate with isst
    void integrate_with_isst(Molecule* molecule_pt, const double& tmin,
                             const double& tmax, const unsigned& nint,
                             const double& t_step, const double& tau);
    // set with grid
    void integrate_with_grid(const double& a_x, const double& b_x,
                             const double& b_y, const double& cut_off,
                             Molecule* molecule_pt);
    // update the pressure and temperature reading
    void npt_update_pressure_temperature();
    double npt_get_instant_pressure();
    // get the pressure reading
    double npt_get_pressure();
    // get the instant temperature reading
    double npt_get_instant_temperature();
    // get the temperature reading
    double npt_get_temperature();
    // set the  initial condition from file
    void npt_set_initial(Molecule* molecule_pt,
                         const char* initial_pos_filename,
                         const char* initial_mom_filename,
                         const char* initial_box_filename);


protected:
    void A(Particle* particle_pt, const double& h);
    void B(Particle* particle_pt, const double& h);
    void O(Particle* particle_pt);

    void A_1_NPT(const double& h);
    void A_2_NPT(Molecule* molecule_pt, const double& h);
    void B_NPT(Molecule* molecule_pt, const double& h);
    void O_NPT(Molecule* molecule_pt);

    // compute the force in the correct way
    void compute_force(Molecule* molecule_pt);
    // set with npt grid
    void integrate_with_npt_grid(const double& a_x, const double& b_x,
                                 const double& b_y, const double& cut_off,
                                 Molecule* molecule_pt, const double& mass,
                                 const double& target_press,
                                 const double& gamma_npt,
                                 const double& o_box_time_step);
    // set to integrate with simulated tempering
    void integrate_with_st(const double& tmin,
                           const double& tmax,
                           const double& n_temperatures,
                           const unsigned& mod_switch,
                           const int& seed);
    // check if to update temperature
    void update_simulated_tempering(Molecule* molecule_pt,
                                    const unsigned& step,
                                    const double& h);

private:
    double Beta;
    double Gamma;
    double GammaRot;
    double GammaNpt;

    double OstepC;
    double OstepZ;
    double OstepCphi;
    double OstepZphi;
    double OstepCbox;
    double OstepZbox;

    double Target_pressure;
    double Npt_mass;

    bool With_isst;
    bool With_grid;
    bool With_npt;

    // gsl_rng* pos_gen;
    // gsl_rng* rot_gen;
    // gsl_rng* box_gen;
    // std::mt19937_64 pos_gen;
    // std::mt19937_64 rot_gen;
    // std::mt19937_64 box_gen;
    NormalGenerator normal_gen;
    //NormalGenerator* rot_gen_pt;
    //NormalGenerator* box_gen_pt;

    // class pointers
    System* System_pt;
    InfiniteSwitchSimulatedTempering* Isst_pt;
    Grid* Grid_pt;
    NptGrid* NptGrid_pt;
    SimulatedTempering* St_pt;

    // rotation steps
    void A_rot(Particle* particle_pt, const double& h);
    void B_rot(Particle* particle_pt, const double& h);
    void O_rot(Particle* particle_pt);

    // NPT helper steps
    void A_NPT_2_part(Particle* particle_pt, const double& h);
    void A_NPT_2_box(const double& h);

    void B_NPT_part(Particle* particle_pt, const double& h);
    void B_NPT_box(const double& h);

    void O_NPT_box();
};

// /******************************************************************************
//                                GENERIC INTEGRATOR
//  *****************************************************************************/
//
// template <class SYSTEM> class Integrator
// {
//  protected:
//
//   // pointer to the system
//   SYSTEM* System_pt;
//
//   // time step
//   double Time_Step;
//
//   // Langevin Friction
//   double Gamma;
//   double Gamma_Rot;
//
//   // Normal distributed random number
//   const gsl_rng_type* T;
//   gsl_rng* r;
//
//   gsl_rng* r_phi;
//
//   // boolean and upper and lower
//   // boundary for the fleming-viot process
//   bool WITH_FLEMING_VIOT;
//
//   bool WITH_TRANS_STATS;
//
//   double fleming_viot_lower;
//   double fleming_viot_upper;
//
//   unsigned fleming_viot_npart;
//   gsl_rng* fleming_viot_index;
//
//   unsigned fleming_viot_tau;
//   unsigned fleming_viot_count;
//
//   Observables* ObservablePt;
//
//   long unsigned TIMEINDEX;
//
//   long unsigned RightTransitionCount;
//   long unsigned LeftTransitionCount;
//
//   // constants for the ostep
//   double Ostep_c;
//   double Ostep_zeta;
//
//   double Ophi_step_c;
//   double Ophi_step_zeta;
//
//   // A step - molecule
//   void A_step(Molecule* Mol_pt, const double& step_size)
//   {
//     // dereference
//     unsigned nParticle = Mol_pt->nparticle();
//     unsigned dim       = Mol_pt->dim();
//
//     double* q;
//     double  p = 0.0;
//     double  m = 0.0;
//
//     // loop over the number of particles
//     for(unsigned i=0;i<nParticle;i++)
//       {
//   	// loop over the number of dimensions
//   	for(unsigned j=0;j<dim;j++)
//   	  {
// 	    q =  Mol_pt->particle_pt(i)->q_pt(j);
// 	    p = *Mol_pt->particle_pt(i)->p_pt(j);
// 	    m = *Mol_pt->particle_pt(i)->m_pt(j);
//
//   	    *q += step_size * Time_Step * p/m;
//   	  }
//       }
//   };
//
//   // B step - molecule
//   void B_step(Molecule* Mol_pt,const double& step_size)
//   {
//     unsigned nParticle = Mol_pt->nparticle();
//     unsigned dim       = Mol_pt->dim();
//
//     double* p;
//     double  f = 0.0;
//
//     // loop over the number of particles
//     for(unsigned i=0;i<nParticle;i++)
//       {
//   	// loop over the number of dimensions
//   	for(unsigned j=0;j<dim;j++)
//   	  {
// 	    p =  Mol_pt->particle_pt(i)->p_pt(j);
// 	    f = *Mol_pt->particle_pt(i)->f_pt(j);
//
// 	    *p += step_size * Time_Step * f;
//   	  }
//       }
//   };
//
//   // O step - molecule
//   void O_step(Molecule* Mol_pt, const double& step_size)
//   {
//     unsigned nParticle = Mol_pt->nparticle();
//     unsigned dim       = Mol_pt->dim();
//
//     double* p;
//     double  m  = 0.0;
//     double  kT = Mol_pt->kt();
//
//     // loop over the number of particles
//     for(unsigned i=0;i<nParticle;i++)
//       {
//   	// loop over the number of dimensions
//   	for(unsigned j=0;j<dim;j++)
//   	  {
// 	    p =  Mol_pt->particle_pt(i)->p_pt(j);
// 	    m = *Mol_pt->particle_pt(i)->m_pt(j);
//
// 	    *p = Ostep_c * (*p) + Ostep_zeta * sqrt(m)
// 	      * gsl_ran_gaussian(r,1.0);
//   	  }
//       }
//   };
//
//   // A step - particle
//   void A_step(Particle* part,
// 	      const double& step_size,
// 	      const unsigned& dim)
//   {
//     double* q;
//     double  p = 0.0;
//     double  m = 0.0;
//
//     // loop over the number of dimensions
//     for(unsigned j=0;j<dim;j++)
//       {
// 	q =  part->q_pt(j);
// 	p = *part->p_pt(j);
// 	m = *part->m_pt(j);
//
// 	*q += step_size * Time_Step * p/m;
//       }
//   };
//
//   // B step - particle
//   void B_step(Particle* part,
// 	      const double& step_size,
// 	      const unsigned& dim)
//   {
//     double* p;
//     double  f = 0.0;
//
//     // loop over the number of dimensions
//     for(unsigned j=0;j<dim;j++)
//       {
// 	p =  part->p_pt(j);
// 	f = *part->f_pt(j);
//
// 	*p += step_size * Time_Step * f;
//       }
//   };
//
//   // O step - particle
//   void O_step(Particle* part,
// 	      const double& step_size,
// 	      const double& kT,
// 	      const unsigned& dim)
//   {
//     double* p;
//     double  m  = 0.0;
//
//     // loop over the number of dimensions
//     for(unsigned j=0;j<dim;j++)
//       {
// 	p =  part->p_pt(j);
// 	m = *part->m_pt(j);
//
// 	*p = Ostep_c * (*p) + Ostep_zeta
// 	  * sqrt(m) *gsl_ran_gaussian(r,1.0);
//       }
//
//   };
//
//   // Bphi step, changes angular momentum
//   void BpPhi_step(Particle* part, const double& step_size,
//                   const unsigned& rot_dim)
//   {
//       double* pi = NULL;
//       double tau = 0.0;
//
//       for(unsigned i=0; i<rot_dim; i++)
//         {
//             pi = part->pi_pt(i);
//             tau = *part->tau_pt(i);
//
//             *pi += step_size * Time_Step * tau;
//         }
//   };
//
//   // ApPhi step, changes angular momentum
//   void ApPhi_step(Particle* part,const double& step_size,
//                   const unsigned& rot_dim)
//   {
//       double pi = 0.0;
//       double I = 0.0;
//       double alpha = 0.0;
//
//       // rotation matrix
//       double R_00 = 0.0;
//       double R_01 = 0.0;
//
//       double R_10 = 0.0;
//       double R_11 = 0.0;
//
//       for(unsigned i=0; i<rot_dim; i++)
//       {
//           pi = *part->pi_pt(i);
//           I = *part->I_pt(i);
//
//           alpha = step_size * Time_Step * pi / I;
//
//           // R = [  Cos(alpha) Sin(alpha)
//           //       -Sin(alpha) Cos(alpha) ]
//           R_00 =  cos(alpha);
//           R_01 = -sin(alpha);
//
//           R_10 =  sin(alpha);
//           R_11 =  cos(alpha);
//
//           // Update the rotation matrix
//           //     Q = QR^T
//           double Q_00 = *part->Q_pt(0,0) * R_00 + *part->Q_pt(0,1) * R_01;
//           double Q_01 = *part->Q_pt(0,0) * R_10 + *part->Q_pt(0,1) * R_11;
//           double Q_10 = *part->Q_pt(1,0) * R_00 + *part->Q_pt(1,1) * R_01;
//           double Q_11 = *part->Q_pt(1,0) * R_10 + *part->Q_pt(1,1) * R_11;
//
//           // Asign the new matrix
//           *part->Q_pt(0,0) = Q_00;
//           *part->Q_pt(0,1) = Q_01;
//
//           *part->Q_pt(1,0) = Q_10;
//           *part->Q_pt(1,1) = Q_11;
//       }
//  };
//
//   // Bphi step, changes angular momentum
//   void OpPhi_step(Particle* part, const double& kT)
//   {
// 	double* pi =  part->pi_pt(0);
//     double I = *part->I_pt(0);
//
//     *pi = Ophi_step_c * (*pi) + Ophi_step_zeta * sqrt(I)
//           * gsl_ran_gaussian(r_phi, 1.0);
//   };
//
//   void initialise_integrator(const double& time_step,
// 			     const double& gamma,
// 			     const double& kT,
// 			     const double& Ostep_size,
//                  const double& gamma_rot)
//   {
//     // Set the time step and Langevin Friction
//     Time_Step = time_step;
//
//     // set the gamma's for rotation and for translation
//     Gamma = gamma;
//     Gamma_Rot = gamma_rot;
//
//     // create pointer to system
//     System_pt = new SYSTEM;
//
//     // set the timeindex to zero
//     TIMEINDEX = 0;
//
//     // Generate a random seed based on system variables
//     long unsigned seed;
//     struct timeval tv;
//     gettimeofday(&tv,0);
//     seed = tv.tv_sec + tv.tv_usec + 15.0 * omp_get_thread_num();
//
//     // Select and seed the random number generator
//     r = gsl_rng_alloc (gsl_rng_mt19937);
//     gsl_rng_set(r, seed);
//
//     gettimeofday(&tv,0);
//     seed = tv.tv_sec + tv.tv_usec + 5.0 * omp_get_thread_num();
//
//     r_phi = gsl_rng_alloc (gsl_rng_mt19937);
//     gsl_rng_set(r_phi, seed);
//
//     // set all the o step parameters
//     update_o_step_parameters(kT, Ostep_size);
//   };
//
//   void update_o_step_parameters(const double& kT, const double& Ostep_size)
//   {
//       // set up the o-step constatns
//       Ostep_c    = exp(-Gamma * Ostep_size * Time_Step);
//       Ostep_zeta = sqrt((1.0 - Ostep_c * Ostep_c) * kT);
//
//       // warning set to BAOAB
//       Ophi_step_c = exp(-Gamma_Rot * Ostep_size * Time_Step);
//       Ophi_step_zeta = sqrt((1.0 - Ophi_step_c * Ophi_step_c) * kT);
//   }
//
//  public:
//   Integrator(){};
//
//   // pure virtual Integrate function
//   virtual void Integrate(Molecule* Mol_pt) = 0;
//
//   // return time step
//   double time_step(){return Time_Step;}
//
//  // return the Langevin Friction
//   double gamma(){return &Gamma;}
//
//   // pointer to the system
//   SYSTEM* system_pt(){return System_pt;}
//
//   // return uniform random variable
//   double uniform(){return gsl_rng_uniform(Integrator<SYSTEM>::r);}
//
//   // set the boundary condition to be Flemming-viot
//   void set_bc_fleming_viot
//   (Molecule* MolPt,
//    const double& lower_boundary,
//    const double& upper_boundary,
//    const unsigned& number_of_particles,
//    const unsigned& tau)
//   {
//     // set the boolean to true
//     WITH_FLEMING_VIOT = true;
//
//     // set the lower and upper bounds
//     fleming_viot_lower = lower_boundary;
//     fleming_viot_upper = upper_boundary;
//
//     // set the total number of particles used in simulation
//     fleming_viot_npart = number_of_particles;
//
//     // set the random number generation
//     long unsigned seed;
//     struct timeval tv;
//     gettimeofday(&tv,0);
//     seed = tv.tv_sec + tv.tv_usec + 15.0 * omp_get_thread_num();
//
//     // Select and seed the random number generator
//     fleming_viot_index = gsl_rng_alloc (gsl_rng_mt19937);
//     gsl_rng_set(fleming_viot_index,seed);
//
//     // set fleming viout count
//     fleming_viot_tau = tau;
//     fleming_viot_count = 0;
//
//     // initalise the interpolation in the molecule
//     // so that we can determine barrier crossings more exactly.
//     // set only two previous point for now.
//     MolPt->initialise_interpolation(2);
//   }
//
//   int detect_fleming_viot(const double& q)
//   {
//     int index = 0;
//
//     // detect if the lower or upper barrier was crossed
//     if(q <= fleming_viot_lower || q > fleming_viot_upper)
//       index = gsl_rng_uniform_int(fleming_viot_index,fleming_viot_npart);
//     else
//       index = -1;
//
//     // return the index
//     return index;
//   };
//
//   // returnw the momentum at the time the barrier was crossed
//   // this can be used to collect statistics with
//   void enforce_fleming_viot(Molecule* MolPt)
//   {
//     // number of particles
//     unsigned nParticle = MolPt->nparticle();
//
// #pragma omp for schedule(static)
//     for(unsigned i=0;i<nParticle;i++)
//       {
// 	// reference particle
// 	Particle* PartPt = MolPt->particle_pt(i);
//
// 	// dereference the position
// 	double* q = PartPt->q_pt(0);
// 	double* p = PartPt->p_pt(0);
// 	int index = 0;
//
// 	// get the momentum at the barrier
// 	// for collecting statistics
// 	double p_star = 0.0;
//
// 	while(index != -1)
// 	  {
// 	    // get the index of the paricle
// 	    index = detect_fleming_viot(*q);
//
// 	    // make a copy of the new particle
// 	    if(index != -1)
// 	      {
// 		// get momentum at the barrier crossing
// 		if(*q <= fleming_viot_lower)
// 		  p_star = PartPt->get_threepnt_interpolation_momentum(*q,fleming_viot_lower,*p,Time_Step);
// 		else if(*q > fleming_viot_upper)
// 		  p_star = PartPt->get_threepnt_interpolation_momentum(*q,fleming_viot_upper,*p,Time_Step);
//
// 		// record the momentum
// 		if(ObservablePt != NULL && TIMEINDEX > ObservablePt->get_burn_in())
// 		  ObservablePt->addsample_momentum_dist(p_star);
//
// 		// copy the postion, force and momentum of the index particle
// 		PartPt->copy_particle(MolPt->particle_pt(index));
//
// 		// copy the history of the index particle
// 		PartPt->copy_interpolation(MolPt->particle_pt(index));
// 	      }
// 	  }
// 	// update the particle interpolation for
// 	// all the particles
// 	MolPt->update_particle_interpolation(i,*q,*p);
//       }
//   };
//
//   void set_observable_pt(Observables* obspt)
//   {
//     // set observable pointer
//     ObservablePt = obspt;
//   };
//
//   void set_trans_stats
//   (Molecule* MolPt,
//    const size_t& N,
//    const double& lower_boundary,
//    const double& upper_boundary,
//    const size_t& Np,
//    const double& p_lower_boundary,
//    const double& p_upper_boundary)
//   {
//     // set the boolean to true
//     WITH_TRANS_STATS = true;
//
//     // set the lower and upper bounds
//     fleming_viot_lower = lower_boundary;
//     fleming_viot_upper = upper_boundary;
//
//     // initialise the transition counters
//     LeftTransitionCount  = 0;
//     RightTransitionCount = 0;
//
//     // initialise the particles local histogram
//     for(unsigned i=0;i<MolPt->nparticle();i++)
//       MolPt->particle_pt(i)->initialise_barrier_statistics(N,
// 							   lower_boundary,
// 							   upper_boundary,
// 							   Np,
// 							   p_lower_boundary,
// 							   p_upper_boundary);
//   }
//
//   void print_transition_counts()
//   {
//     printf("\\ ----------------------------------------- \\ \n");
//     printf("Number of Transition Left  to Right: %lu\n",LeftTransitionCount);
//     printf("Number of Transition Right to Left: %lu\n",RightTransitionCount);
//     printf("\\ ----------------------------------------- \\ \n");
//   }
//
//   void collect_trans_stat(Molecule* MolPt)
//   {
//     // number of particles
//     unsigned nParticle = MolPt->nparticle();
//
//     // loop over the particles
// #pragma omp for schedule(static)
//     for(unsigned i=0;i<nParticle;i++)
//       {
// 	// reference particle
// 	Particle* PartPt = MolPt->particle_pt(i);
//
// 	// dereference the position
// 	double q = *PartPt->q_pt(0);
//
// 	// STAGE 1 - detect and update the region
//
// 	// left region
// 	if(q < fleming_viot_lower && PartPt->get_label(0) != -1)
// 	  {
// 	    // set the label of particle
// 	    PartPt->set_label(0,-1);
//
// 	  } // trans region
// 	else if(fleming_viot_lower <= q && q < fleming_viot_upper)
// 	  {
// 	    // set the label of particle unless
// 	    // its already zero
// 	    if(PartPt->get_label(0) != 0)
// 	      PartPt->set_label(0,0);
//
// 	    // add to particle statistics
// 	    PartPt->addsample_barrier_statistics();
//
// 	  } // upper region
// 	else if(fleming_viot_upper < q && PartPt->get_label(0) != 1)
// 	  {
// 	    // set the label of particle
// 	    PartPt->set_label(0,1);
// 	  }
//
// 	// dereference the labels of the particle
// 	int LabelZero = PartPt->get_label(0);
// 	int LabelOne  = PartPt->get_label(1);
//
// 	// STAGE 2 - based on region we found decide what to do
//
// 	// Transiton to region 1 from region -1
// 	if(LabelZero == 1 && LabelOne == -1)
// 	  {
// 	    // dump the statistics into histogram
// 	    if(ObservablePt != NULL && TIMEINDEX > ObservablePt->get_burn_in())
// 	      {
// 		ObservablePt->addsample_pdf_dist(PartPt->get_barrier_q_statistics());
// 		ObservablePt->addsample_momentum_dist(PartPt->get_barrier_p_statistics());
// 	      }
//
// 	    // add to the transition variable
// #pragma omp atomic
// 	    RightTransitionCount++;
//
// 	    // reset the statistics on the particle
// 	    PartPt->reset_barrier_statistics();
//
// 	    // set previos label of the particle to 1
// 	    PartPt->set_label(1,1);
//
// 	  } // Transition from region 1 to -1
// 	else if(LabelZero == -1 && LabelOne == 1)
// 	  {
// 	    // dump statistics into histogram
// 	    if(ObservablePt != NULL && TIMEINDEX > ObservablePt->get_burn_in())
// 	      {
// 		ObservablePt->addsample_pdf_dist(PartPt->get_barrier_q_statistics());
// 		ObservablePt->addsample_momentum_dist(PartPt->get_barrier_p_statistics());
// 	      }
//
// 	    // add to the transition variable
// #pragma omp atomic
// 	    LeftTransitionCount++;
//
// 	    // reset the statistics on the particle
// 	    PartPt->reset_barrier_statistics();
//
// 	    // set the previos label of the particle to -1
// 	    PartPt->set_label(1,-1);
//
// 	  } // No Transition left region 1 and came back to region 1
// 	else if(LabelZero == 1 && LabelOne == 1)
// 	  {
// 	    // reset the statistics on the particle
// 	    PartPt->reset_barrier_statistics();
//
// 	  } // No Transition left region -1 and came back to region -1
// 	else if(LabelZero == -1 && LabelOne == -1)
// 	  {
// 	    // reset the statistics on the particle
// 	    PartPt->reset_barrier_statistics();
// 	  }
//       }
//   };
//
//   void set_initial_condition_exact(Molecule* MolPt,
// 				   const double& a,
// 				   const double& b)
//   {
//     // number of particles
//     unsigned nParticle = MolPt->nparticle();
//     double beta = 1.0/MolPt->kt();
//
//     // set the random number generation
//     long unsigned seed;
//     struct timeval tv;
//     gettimeofday(&tv,0);
//     seed = tv.tv_sec + tv.tv_usec + 15.0 * omp_get_thread_num();
//
//     // Select and seed the random number generator
//     fleming_viot_index = gsl_rng_alloc (gsl_rng_mt19937);
//     gsl_rng_set(fleming_viot_index,seed);
//
//     Particle* PartPt;
//     double alpha = -1.0;
//     double prop = 0.0;
//     double u = 0.0;
//
//     // loop over the particles
// #pragma omp parallel for schedule(static) firstprivate(alpha,prop,u) private(PartPt)
//     for(unsigned i=0;i<nParticle;i++)
//       {
// 	// reference particle
// 	PartPt = MolPt->particle_pt(i);
//
// 	// calculate acceptance
// 	alpha = -1.0;
//
// 	prop = 0.0;
// 	u = 0.0;
//
// 	while(alpha < 0.0)
// 	  {
// 	    // generate proposal
// 	    prop = gsl_ran_flat(fleming_viot_index,a,b);
//
// 	    // calculate acceptance
// 	    alpha = System_pt->get_exact_pdf(prop,1.0/MolPt->kt())/ 3.0;
//
// 	    // generate uniform
// 	    u = gsl_ran_flat(fleming_viot_index,0.0,1.0);
//
// 	    // accept or reject
// 	    if(u <= alpha)
// 	      {
// 		*PartPt->q_pt(0) = prop;
// 		*PartPt->p_pt(0) = sqrt(1.0/beta) * gsl_ran_gaussian(r,1.0);
// 	      }
// 	    else
// 	      alpha = -1.0;
// 	  }
//
//       }
//   };
//
//
//   void set_pos_random(Molecule* MolPt,
// 		      const double& mean,
// 		      const double& var)
//   {
//     // get paramaters
//     unsigned N   = MolPt->nparticle();
//     unsigned DIM = MolPt->dim();
//
//     // loop over all the particles
//     for(unsigned i=0;i<N;i++)
//       for(unsigned j=0;j<DIM;j++)
// 	*MolPt->particle_pt(i)->q_pt(j) = mean + gsl_ran_gaussian(r,sqrt(var));
//   };
//
//   void set_temperature(Molecule* MolPt)
//   {
//     // get paramaters
//     unsigned N   = MolPt->nparticle();
//     unsigned DIM = MolPt->dim();
//
//     double beta = 1.0/MolPt->kt();
//
//     // loop over all the particles
//     for(unsigned i=0;i<N;i++)
//       for(unsigned j=0;j<DIM;j++)
// 	*MolPt->particle_pt(i)->p_pt(j) = sqrt(1.0/beta) * gsl_ran_gaussian(r,1.0);
//   };
// };
//
// /******************************************************************************
//                         BAOABpfc INTEGRATOR - parallel force cell
//  *****************************************************************************/
// template <class SYSTEM> class BAOABpfc : public Integrator<SYSTEM>
// {
//  private:
//
//   // pointer to the CellGrid
//   Grid* CellGridPt;
//
//   // The Total Box is [-Lx, Lx, -Ly, Ly]
//   double Lx;
//   double Ly;
//
//   // Holder for the cutoff in the force calculation
//   double CutOff;
//
//   // Holder for the number of boxes
//   unsigned TotNumBoxes;
//
//  public:
//   // Constructor
//   BAOABpfc(const double& time_step,
// 	   const double& gamma)
//     {
//       // initialise the integrator
//       Integrator<SYSTEM>::initialise_integrator(time_step,gamma);
//     };
//
//   void pre_force_integration(Particle* PartPt,
// 			     const unsigned& dim,
// 			     const double& kT)
//   {
//     Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
//     Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
//     Integrator<SYSTEM>::O_step(PartPt,1.0,dim,kT);
//     Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
//
//     if(dim != 1)
//       {
// 	// Integrate angular momentum forward
// 	Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
//
// 	// Integrate rotation forward
// 	Integrator<SYSTEM>::ApPhi_step(PartPt);
//       }
//   };
//
//   void post_force_integration(Particle* PartPt,
// 			      const unsigned& dim)
//   {
//     Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
//
//     if(dim != 1)
//       {
// 	// integrate angular momentum forward
// 	Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
//       }
//   };
//
//   // Integrator function
//   void Integrate(Molecule* Mol_pt)
//   {
//  /*    // Holder for calculating the potential for each timestep */
// /*     double V = 0.0; */
//
// /* #pragma omp parallel */
// /*     { */
//
// /* #pragma omp for schedule(static) */
// /*       for(unsigned i=0;i<*Mol_pt->nparticle_pt();i++) */
// /* 	{ */
// /* 	  Bp_step(Mol_pt->particle_pt(i),0.5); */
// /* 	  Ap_step(Mol_pt->particle_pt(i),0.5); */
// /* 	  //Op_step(Mol_pt->particle_pt(i),1.0,Mol_pt->kt()); */
// /* 	  Ap_step(Mol_pt->particle_pt(i),0.5); */
//
// /* 	  // Integrate angular momentum forward */
// /* 	  BpPhi_step(Mol_pt->particle_pt(i),0.5); */
//
// /* 	  // Integrate rotation forward */
// /* 	  ApPhi_step(Mol_pt->particle_pt(i)); */
// /* 	} */
//
//
// /*     // update the positions of the particles in the boxes */
// /* #pragma omp single */
// /*       { */
// /* 	Integrator<SYSTEM>::System_pt->grid_pt()->update_particles_on_grid(Mol_pt); */
// /* 	*Mol_pt->potential_pt()=0.0;	 */
// /*       } */
// /* #pragma omp barrier */
//
// /*       // Force Solve - loop over boxes */
// /* #pragma omp for schedule(static) reduction(+:V) */
// /*       for(unsigned i=0;i<Integrator<SYSTEM>::System_pt->get_total_box_num();i++) */
// /* 	{ */
// /* 	  Integrator<SYSTEM>::System_pt->compute_force(Mol_pt,i,V); */
// /* 	} */
//
// /* #pragma omp barrier */
// /* #pragma omp single */
// /*       { */
// /* 	*Mol_pt->potential_pt() = V;	 */
// /*       } */
//
// /*       // final B step for new force level */
// /* #pragma omp for schedule(static) */
// /*       for(unsigned i=0;i<*Mol_pt->nparticle_pt();i++) */
// /* 	{ */
// /* 	  Bp_step(Mol_pt->particle_pt(i),0.5); */
//
// /* 	  // integrate angular momentum forward */
// /* 	  BpPhi_step(Mol_pt->particle_pt(i),0.5); */
//
// /* 	} */
// /*     } */
//   }
//
//   // Destructor
//   ~BAOABpfc()=delete;
// };
//
// /******************************************************************************
//                                 BAOAB INTEGRATOR
//  *****************************************************************************/
// // defined in seperate file.
//
// /******************************************************************************
//                                 ABOBA INTEGRATOR
//  *****************************************************************************/
// // template <class SYSTEM> class ABOBA : public Integrator<SYSTEM>
// // {
// //  public:
// //   // Constructor
// //   ABOBA(const double& time_step,const double& gamma)
// //     {
// //       // Set the time step and Langevin Friction
// //       Integrator<SYSTEM>::Time_Step=time_step;
// //       Integrator<SYSTEM>::Gamma=gamma;
//
// //       // create pointer to system
// //       Integrator<SYSTEM>::System_pt=new SYSTEM;
//
// //       // Generate a random seed based on system variables
// //       long unsigned seed;
// //       struct timeval tv;
// //       gettimeofday(&tv,0);
// //       seed = tv.tv_sec + tv.tv_usec;
//
// //       // Select and seed the random number generator
// //       Integrator<SYSTEM>::r = gsl_rng_alloc (gsl_rng_mt19937);
// //       gsl_rng_set(Integrator<SYSTEM>::r,seed);
// //     };
//
// //   // Integrator function
// //   void Integrate(Molecule* Mol_pt)
// //   {
// //     // F is force from previous step
// //     // The first steps at this force level
// //     Integrator<SYSTEM>::A_step(Mol_pt,0.5);
//
// //     // Force Solve
// //     Integrator<SYSTEM>::System_pt->compute_force(Mol_pt);
//
// //     // final steps for new force level
// //     Integrator<SYSTEM>::B_step(Mol_pt,0.5);
// //     Integrator<SYSTEM>::O_step(Mol_pt,1.0);
// //     Integrator<SYSTEM>::B_step(Mol_pt,0.5);
// //     Integrator<SYSTEM>::A_step(Mol_pt,0.5);
// //   }
// //   // Destructor
// //   ~ABOBA()=delete;
// // };
//
// /******************************************************************************
//                                 OBABO INTEGRATOR
//  *****************************************************************************/
// template <class SYSTEM> class OBABO : public Integrator<SYSTEM>
// {
//  public:
//   // Constructor
//   OBABO(const double& time_step,const double& gamma)
//     {
//       // Set the time step and Langevin Friction
//       Integrator<SYSTEM>::Time_Step=time_step;
//       Integrator<SYSTEM>::Gamma=gamma;
//
//       // create pointer to system
//       Integrator<SYSTEM>::System_pt=new SYSTEM;
//
//       // Generate a random seed based on system variables
//       long unsigned seed;
//       struct timeval tv;
//       gettimeofday(&tv,0);
//       seed = tv.tv_sec + tv.tv_usec;
//
//       // Select and seed the random number generator
//       Integrator<SYSTEM>::r = gsl_rng_alloc (gsl_rng_mt19937);
//       gsl_rng_set(Integrator<SYSTEM>::r,seed);
//     };
//
//   // Integrator function
//   void Integrate(Molecule* Mol_pt)
//   {
//     // F is force from previous step
//     // The first steps at this force level
//     Integrator<SYSTEM>::O_step(Mol_pt,0.5);
//     Integrator<SYSTEM>::B_step(Mol_pt,0.5);
//     Integrator<SYSTEM>::A_step(Mol_pt,1.0);
//
//     // Force Solve
//     Integrator<SYSTEM>::System_pt->compute_force(Mol_pt);
//
//     // final steps for new force level
//     Integrator<SYSTEM>::B_step(Mol_pt,0.5);
//     Integrator<SYSTEM>::O_step(Mol_pt,0.5);
//   }
//   // Destructor
//   ~OBABO()=delete;
// };
//
// /* /\****************************************************************************** */
// /*                            IREMD_BAOAB INTEGRATOR */
// /*  *****************************************************************************\/ */
// /* template <class SYSTEM> class IREMD_BAOAB : public Integrator<SYSTEM> */
// /* { */
// /*  private: */
//
// /*   // Infinite REMD B step, changes momentum */
// /*   void IREMD_B_step(Molecule* Mol_pt,const double& step_size) */
// /*   { */
// /*     // loop over the number of particles */
// /*     for(unsigned i=0;i<*Mol_pt->nparticle_pt();i++) */
// /*       { */
// /*   	// loop over the number of dimensions */
// /*   	for(unsigned j=0;j<*Mol_pt->dim_pt();j++) */
// /*   	  { */
// /*   	    *Mol_pt->particle_pt(i)->p_pt(j)=*Mol_pt->particle_pt(i)->p_pt(j) */
// /* 	      +step_size*(Integrator<SYSTEM>::Time_Step) */
// /* 	      *(*Integrator<SYSTEM>::System_pt->r_pt())* */
// /* 	      (*Mol_pt->particle_pt(i)->f_pt(j)); */
// /*   	  } */
// /*       } */
// /*   }; */
//
// /*   // Calculate the Rk's */
// /*   void Calculate_scaling(std::vector<double>& potentials, */
// /* 			 std::vector<double>& R,  */
// /* 			 std::vector<double>& beta, */
// /* 			 std::vector<double>& Omega) */
// /*   { */
// /*     // temporary doubles */
// /*     double tmp=1.0; */
// /*     double sumk =0.0; */
// /*     unsigned i=0; */
//
// /*     // determine the number of temperatures */
// /*     unsigned N=potentials.size();  */
//
// /*     // Loop to determine the sum of permutations for R_k */
// /*     for (unsigned k=0; k<N;k++) */
// /*       { */
// /* 	// Clear the previous force scaling */
// /* 	R[k]=0.0; */
//
// /* 	// loop for all the permutations */
// /* 	do */
// /* 	  { */
//
// /* 	    // loop over all the current permutations */
// /* 	    for(unsigned j=0; j<N;j++) */
// /* 	      { */
// /* 	    	tmp *= exp(-beta[j]*potentials[j]); */
// /* 	      } */
//
// /* 	    // add the permutations together and sum the tmp variable */
// /* 	    R[k]+=beta[k]*tmp; */
// /* 	    sumk+=tmp; */
//
// /* 	    // save the omegas only on first run. */
// /* 	    if(k==0) */
// /* 	      { */
// /* 	    	Omega[i]=tmp; */
// /* 	    	i++; */
// /* 	      } */
//
// /* 	    // reset the tmp variable */
// /* 	    tmp=1.0; */
//
// /* 	  } while(std::next_permutation(beta.begin(),beta.end())); */
//
// /* 	// divide R_k and reset the sumk */
// /* 	R[k]/=(sumk*beta.back()); */
//
// /* 	// divide by the sum of all to get the omegas,this only works if we */
// /* 	// have 2 temperatures */
// /* 	if(k==0) */
// /* 	  { */
// /* 	    for(unsigned l=0;l<Omega.size();l++) */
// /* 	      { */
// /* 		Omega[l]/=sumk; */
// /* 	      } */
// /* 	  } */
// /* 	// reset sum */
// /* 	sumk = 0.0; */
// /*       } */
// /*   }; */
//
// /*   // Parallel step in integrator - calculates either Force or the omegas */
// /*   // depending on the rank of the processor */
// /*   void Parallel_step(Molecule* Mol_pt,double omega[], */
// /* 		     std::vector<double>& beta, */
// /* 		     double potential[], double R[], */
// /* 		     const unsigned& N, const unsigned& id) */
// /*     { */
// /*       if(id==0) */
// /* 	{ */
// /* 	  // first remove the first entry which is dummy */
// /* 	  std::vector<double> pot_tmp(N,0.0); */
// /* 	  std::vector<double> R_tmp(N,0.0); */
// /* 	  std::vector<double> Omega_tmp(factorial(N),0.0); */
//
// /* 	  for(unsigned i=1;i<N+1;i++) */
// /* 	    { */
// /* 	      pot_tmp[i-1]=potential[i]; */
// /* 	    } */
//
// /* 	  // calculate the R's */
// /* 	  Calculate_scaling(pot_tmp,R_tmp,beta,Omega_tmp); */
//
// /* 	  // set dummy R[0] to zero */
// /* 	  R[0]=0.0; */
// /* 	  omega[0]=0.0; */
//
// /* 	  // update the correct R's  */
// /* 	  for(unsigned i=1;i<N+1;i++) */
// /* 	    { */
// /* 	      R[i]=R_tmp[i-1]; */
// /* 	    } */
//
// /* 	  // update the omegas */
// /* 	  for(unsigned i=1;i<=factorial(N);i++) */
// /* 	    { */
// /* 	      omega[i]=Omega_tmp[i-1]; */
// /* 	    } */
// /* 	} */
// /*       else */
// /*       	{ */
// /*       	  // Calculate the force for all the other  */
// /*       	  Integrator<SYSTEM>::System_pt->compute_force(Mol_pt); */
// /*       	} */
// /*     }; */
//
//
// /*  public:  */
// /*   // Constructor */
// /*   IREMD_BAOAB(const double& time_step,const double& gamma) */
// /*     { */
// /*       // Set the time step and Langevin Friction */
// /*       Integrator<SYSTEM>::Time_Step=time_step; */
// /*       Integrator<SYSTEM>::Gamma=gamma; */
//
// /*       // create pointer to system */
// /*       Integrator<SYSTEM>::System_pt=new SYSTEM; */
//
// /*       // Generate a random seed based on system variables */
// /*       long unsigned seed; */
// /*       struct timeval tv; */
// /*       gettimeofday(&tv,0); */
// /*       seed = tv.tv_sec + tv.tv_usec; */
//
// /*       // Select and seed the random number generator  */
// /*       Integrator<SYSTEM>::r = gsl_rng_alloc (gsl_rng_mt19937); */
// /*       gsl_rng_set(Integrator<SYSTEM>::r,seed); */
// /*     }; */
//
// /*   // Integrator function */
// /*   void Integrate(Molecule* Mol_pt,double omega[], */
// /* 		 std::vector<double>& beta, */
// /* 		 const unsigned& N,const unsigned& id) */
// /*   { */
// /*     // make holders for the potential's and the scaling parameters R */
// /*     double potential[N+1]; */
// /*     double R[N+1]; */
//
// /*     // F is force from previous step. Only integrate if rank isn't zero. */
// /*     if(id!=0) */
// /*       { */
// /* 	IREMD_B_step(Mol_pt,0.5); */
// /* 	Integrator<SYSTEM>::A_step(Mol_pt,0.5); */
// /* 	Integrator<SYSTEM>::O_step(Mol_pt,1.0); */
// /* 	Integrator<SYSTEM>::A_step(Mol_pt,0.5); */
//
// /* 	// Calculate the potential */
// /* 	Integrator<SYSTEM>::System_pt->compute_potential(Mol_pt); */
// /*       } */
//
// /*     // Gather the potentials on process 0 */
// /*     MPI_Gather(Mol_pt->potential_pt(),1,MPI_DOUBLE, */
// /* 	       potential,1,MPI_DOUBLE,0,MPI_COMM_WORLD); */
//
// /*     // perform the parallel step */
// /*     Parallel_step(Mol_pt,omega,beta,potential,R,N,id); */
//
// /*     // scatter the R's around all the replicas */
// /*     MPI_Scatter(R,1,MPI_DOUBLE, */
// /* 		Integrator<SYSTEM>::System_pt->r_pt(), */
// /* 		1,MPI_DOUBLE,0,MPI_COMM_WORLD); */
//
// /*     // integrate last step only if rank isn't 0 */
// /*     if(id!=0) */
// /*       { */
// /* 	// last integration step */
// /* 	IREMD_B_step(Mol_pt,0.5); */
// /*       } */
// /*   } */
//
// /*   // Integrator pre-setup */
// /*   void Parallel_Initialise(Molecule* Mol_pt,double omega[], */
// /* 			   std::vector<double>& beta, */
// /* 			   const unsigned& N,const unsigned& id) */
// /*   { */
// /*     // make holders for the potential's and the scaling parameters R */
// /*     double potential[N+1]; */
// /*     double R[N+1]; */
//
// /*     // Gather all the potentials on processor 0 */
// /*     MPI_Gather(Mol_pt->potential_pt(),1,MPI_DOUBLE, */
// /*     	       potential,1,MPI_DOUBLE,0,MPI_COMM_WORLD); */
//
// /*     if(id==0) */
// /*       { */
// /* 	Parallel_step(Mol_pt,omega,beta,potential,R,N,id); */
// /*       } */
// /*     else */
// /*       { */
// /* 	// do nothing for these id's */
// /*       } */
//
// /*     // scatter the R's around all the replicas */
// /*     MPI_Scatter(R,1,MPI_DOUBLE, */
// /* 		Integrator<SYSTEM>::System_pt->r_pt(), */
// /* 		1,MPI_DOUBLE,0,MPI_COMM_WORLD); */
// /*   }; */
//
//
// /*   // Destructor */
// /*   ~IREMD_BAOAB()=delete; */
// /* }; */
//
// /* /\****************************************************************************** */
// /*                            ISST_BAOAB INTEGRATOR */
// /*  *****************************************************************************\/ */
//
// // Beta integral structure
// struct BetaParams {double potential; double dimension;};
//
// /* // function for calculating beta bar in ISST -- only Harmonic Oscillator */
// /* double fBetaBarNumerator(double beta, void* params); */
//
// /* // function for calculating beta bar in ISST -- only Harmonic Oscillator */
// /* double fBetaBarDeNumerator(double beta, void* params); */
//
// /* template <class SYSTEM> class ISST_BAOAB_OLD : public Integrator<SYSTEM> */
// /* { */
// /*  private: */
// /*   // holder for Physical temperature */
// /*   std::vector<double> Beta; */
// /*   std::vector<double> GaussW; */
//
// /*   // Boolean for if learning should be done or not */
// /*   bool WithLearn; */
// /*   bool OnTheFly; */
//
// /*   // Beta limits for Integration over Beta in BetaBar equation */
// /*   double BetaOne; */
// /*   double BetaTwo; */
//
// /*   // Holder for Alpha and Dimension */
// /*   double Alpha; */
// /*   unsigned Dim; */
//
// /*   // Calculate the termal scaling factor */
// /*   double ThermalScalingFactor(Molecule* Mol_pt) */
// /*   { */
// /*     //initialise factor to return */
// /*     double ThermalScaling = 0.0; */
// /*     unsigned nbeta        = Beta.size(); */
//
// /*     double V              = *Mol_pt->potential_pt(); */
// /*     double BarNumSum      = 0.0; */
//
// /*     double BarDeNumSum    = 0.0; */
//
// /*     std::vector<double> BetaWeight = *Mol_pt->get_beta_weight_pt(); */
//
// /*     //---------------- Update BetaBar -------------------// */
// /*     for(unsigned i=0;i<nbeta;i++) */
// /*       { */
// /* 	BarNumSum   += GaussW[i] * Beta[i] * BetaWeight[i] * exp(-Beta[i] * V); */
// /* 	BarDeNumSum += GaussW[i] * BetaWeight[i] * exp(-Beta[i] * V); */
// /*       } */
//
// /*     ThermalScaling = BarNumSum/((1/Mol_pt->kt()) * BarDeNumSum); */
//
// /*     //--------------- Update BetaWeight -------------------// */
// /*     if(OnTheFly) */
// /*       { */
// /* 	double G_inv = 0.0; */
//
// /* 	for(unsigned i=0;i<nbeta;i++) */
// /* 	  { */
// /* 	    (*Mol_pt->get_R_T_pt())[i] += exp(-Beta[i] * V) * BetaWeight[i] / BarDeNumSum; */
// /* 	    G_inv += GaussW[i] * BetaWeight[i] / ((*Mol_pt->get_R_T_pt())[i]); */
// /* 	  } */
//
// /* 	// update the BetaWeight */
// /* 	for(unsigned n=0;n<nbeta;n++) */
// /* 	  { */
// /* 	    (*Mol_pt->get_beta_weight_pt())[n] = (1.0 - 1.0/Alpha) * BetaWeight[n] */
// /* 	      + 1.0/(Alpha * ((1.0 / BetaWeight[n] * (*Mol_pt->get_R_T_pt())[n]) * G_inv )); */
// /* 	  } */
// /*       } */
// /*     else // !OnTheFly - Not on the fly onlt update R_T vector and leave Weight Alone */
// /*       { */
// /* 	for(unsigned i=0;i<nbeta;i++) */
// /* 	  { */
// /* 	    (*Mol_pt->get_R_T_pt())[i] += exp(-Beta[i] * V) * BetaWeight[i] / BarDeNumSum; */
// /* 	  } */
// /*       } */
//
// /*     return ThermalScaling; */
// /*   } */
//
// /*  public:  */
// /*   // Constructor */
// /*   ISST_BAOAB_OLD(const double& time_step, */
// /* 		 const double& gamma, */
// /* 		 const unsigned& d, */
// /* 		 const double& tmin, */
// /* 		 const double& tmax, */
// /* 		 const double& alpha, */
// /* 		 const unsigned& nInt, */
// /* 		 const bool withlearn, */
// /* 		 const bool onthefly) */
// /*     { */
// /*       // Set the time step and Langevin Friction */
// /*       Integrator<SYSTEM>::Time_Step=time_step; */
// /*       Integrator<SYSTEM>::Gamma=gamma; */
//
// /*       // create pointer to system */
// /*       Integrator<SYSTEM>::System_pt=new SYSTEM; */
//
// /*       // set the boolean variable if with learning */
// /*       WithLearn = withlearn; */
// /*       OnTheFly  = onthefly; */
//
// /*       BetaTwo   = 1.0/tmin; */
// /*       BetaOne   = 1.0/tmax; */
//
// /*       Dim = d; */
//
// /*       if(WithLearn) */
// /* 	{ */
// /* 	  // make the integration table */
// /* 	  gsl_integration_glfixed_table * tbl = gsl_integration_glfixed_table_alloc(nInt); */
//
// /* 	  // dereference the integration points */
// /* 	  for(unsigned i=0;i<nInt;i++) */
// /* 	    { */
// /* 	      // push back to add integration and weight points */
// /* 	      Beta.push_back(0.0); */
// /* 	      GaussW.push_back(0.0); */
//
// /* 	      gsl_integration_glfixed_point(BetaOne,BetaTwo,i,&Beta[i],&GaussW[i],tbl); */
// /* 	    } */
//
// /* 	  // free the integration table */
// /* 	  gsl_integration_glfixed_table_free(tbl); */
//
// /* 	  // initialise Alpha */
// /* 	  Alpha       = alpha; */
// /* 	} */
//
// /*       // Generate a random seed based on system variables */
// /*       long unsigned seed; */
// /*       struct timeval tv; */
// /*       gettimeofday(&tv,0); */
// /*       seed = tv.tv_sec + tv.tv_usec; */
//
// /*       // Select and seed the random number generator  */
// /*       Integrator<SYSTEM>::r = gsl_rng_alloc (gsl_rng_mt19937); */
// /*       gsl_rng_set(Integrator<SYSTEM>::r,seed); */
// /*     }; */
//
// /*   void print_temps(Molecule* MolPt,FILE* output) */
// /*   { */
//
// /*     for(unsigned i=0;i<Beta.size();i++) */
// /*       { */
// /* 	fprintf(output,"%f\t%.10f\n",Beta[i],(*MolPt->get_beta_weight_pt())[i]); */
// /*       }   */
// /*   } */
//
// /*   void print_error(Molecule* MolPt,const unsigned& T, FILE* output) */
// /*   { */
// /*     double error = 0.0; */
// /*     unsigned ntemps = Beta.size(); */
//
// /*     unsigned d = *MolPt->dim_pt(); */
//
// /*     for(unsigned i=0;i<ntemps;i++) */
// /*       { */
// /* 	error += pow(((*MolPt->get_beta_weight_pt())[i] - ((double)d/2.0 + 1) * pow(Beta[i],((double)d/2.0)) */
// /* 		      /(pow(BetaTwo,((double)d/2.0 + 1)) - pow(BetaOne,((double)d/2.0 + 1)))),2.0); */
// /*       } */
//
// /*     // output the abs of the error */
// /*     fprintf(output,"%d\t%.8f\n",T,sqrt(error)); */
// /*   } */
//
// /*   void initialise_beta_weight(Molecule* MolPt,const unsigned& Ic) */
// /*   { */
// /*     unsigned ntemps = Beta.size(); */
//
// /*     double sumbeta = 0.0; */
//
// /*     for(unsigned i=0;i<ntemps;i++) */
// /*       { */
// /* 	sumbeta += Beta[i]; */
// /*       } */
//
// /*     unsigned d = *MolPt->dim_pt(); */
//
// /*     // assign the temperatures */
// /*     for(unsigned i=0;i<ntemps;i++) */
// /*       { */
// /* 	// initialise the R_T vector too -- if this isn't done it produces nan */
// /* 	(*MolPt->get_R_T_pt())[i] = 0.0; */
//
// /* 	switch(Ic) */
// /* 	  { */
// /* 	  case 0: 	 */
// /* 	    // This is the Theoretical value:                   SAVED in BetaTest-exact */
// /* 	    // */
// /* 	    //       Z^-1(b) = b^(d/2) */
// /* 	    // */
// /* 	    (*MolPt->get_beta_weight_pt())[i] = ((double)d/2.0 + 1) * pow(Beta[i],((double)d/2.0)) */
// /* 	      /(pow(BetaTwo,((double)d/2.0 + 1)) - pow(BetaOne,((double)d/2.0 + 1))); */
// /* 	    break; */
//
// /* 	  case 1: */
// /* 	    // set initial condition to power of beta:          SAVED in BetaTest-pow2 */
// /* 	    (*MolPt->get_beta_weight_pt())[i] = pow(10,-6.0)* pow(Beta[i],2.0); */
// /* 	    break; */
//
// /* 	  case 2: */
// /* 	    // set initial conditon to one:                    SAVED in BetaTest-one */
// /* 	    (*MolPt->get_beta_weight_pt())[i] = 1.0; */
// /* 	    break; */
//
// /* 	  case 3: */
// /* 	    // set initial condition to beta:                   SAVED in BetaTest-beta */
// /* 	    (*MolPt->get_beta_weight_pt())[i] = Beta[i]/sumbeta; */
// /* 	    break; */
// /* 	  } */
// /*       } */
// /*   } */
//
// /*   // Integrator function */
// /*   void Integrate(Molecule* Mol_pt) */
// /*   { */
// /*     // F is force from previous step */
// /*     // The first steps at this force level */
// /*     Integrator<SYSTEM>::B_step(Mol_pt,0.5); */
// /*     Integrator<SYSTEM>::A_step(Mol_pt,0.5); */
// /*     Integrator<SYSTEM>::O_step(Mol_pt,1.0); */
// /*     Integrator<SYSTEM>::A_step(Mol_pt,0.5); */
//
// /*     // Force Solve */
// /*     Integrator<SYSTEM>::System_pt->compute_force(Mol_pt); */
// /*     Integrator<SYSTEM>::System_pt->compute_potential(Mol_pt); */
//
// /*     // Evaluate ThermalScaling; */
// /*     double scaling   = 0.0; */
// /*     if(WithLearn) */
// /*       { */
// /* 	scaling    = ThermalScalingFactor(Mol_pt); */
// /*       } */
// /*     else */
// /*       { */
// /* 	scaling    = Integrator<SYSTEM>::System_pt-> */
// /* 	  get_beta_bar(*Mol_pt->potential_pt(), */
// /* 		       *Mol_pt->dim_pt(), */
// /* 		       BetaOne, */
// /* 		       BetaTwo); */
// /*       } */
//
// /*     // rescale the Force */
// /*     for(unsigned i=0;i<*Mol_pt->nparticle_pt();i++) */
// /*       { */
// /* 	for(unsigned j=0;j<*Mol_pt->dim_pt();j++) */
// /* 	  { */
//
// /* 	    *Mol_pt->particle_pt(i)->f_pt(j) *= scaling; */
// /* 	  } */
// /*       } */
//
// /*     // final B step for new force level */
// /*     Integrator<SYSTEM>::B_step(Mol_pt,0.5); */
// /*   } */
// /*   // Destructor */
// /*   ~ISST_BAOAB_OLD()=delete; */
//
// /*   // update the weight functions */
// /*   void update_omega(Molecule* Mol_pt, const unsigned N) */
// /*   { */
// /*     unsigned nbeta = Beta.size(); */
// /*     double G_inv = 0.0; */
//
// /*     // Do G^-1 integral */
// /*     for(unsigned i=0;i<nbeta;i++) */
// /*       { */
// /* 	G_inv += GaussW[i] * (*Mol_pt->get_beta_weight_pt())[i] / ((*Mol_pt->get_R_T_pt())[i]); */
// /*       } */
//
// /*     for(unsigned i=0;i<nbeta;i++) */
// /*       { */
// /* 	(*Mol_pt->get_beta_weight_pt())[i] = ((*Mol_pt->get_R_T_pt())[i]) /  */
// /* 	  (*Mol_pt->get_beta_weight_pt())[i] * G_inv; */
// /*       } */
// /*   } */
//
// /*   // function for returning the number of temperatures */
// /*   unsigned number_of_temps(){return Beta.size();} */
//
// /*   // rewieghting function that returns that adds the correct weight to the correct bin */
// /*   // see paper with Jianfeng, Eric and Ben. Always at lowest temperature */
// /*   void weighted_1d_histogram(const double& A, Molecule* MolPt, gsl_histogram* hist) */
// /*   { */
// /*     // A is the observable which we are after, should be double */
// /*     // R_T is held by the integrator itself and gets initiatied in the constructor */
//
// /*     // need to determine qunatity: */
// /*     // */
// /*     // exp(-beta * V(x)) * omega(beta)/  */
// /*     // int_{beta_1}^{beta_2} exp(-beta * V(x)) * omega(beta) dbeta */
// /*     // */
//
// /*     if(WithLearn) */
// /*       { */
// /* 	// senare */
// /*       } */
// /*     else */
// /*       { */
// /* 	// make holder for the intermideiate weight */
// /* 	double w = Integrator<SYSTEM>::System_pt->get_weight_func(*MolPt->potential_pt(), */
// /* 								  *MolPt->dim_pt(), */
// /* 								  BetaOne, BetaTwo); */
//
// /* 	// add to the holder functions */
// /* 	(*MolPt->get_R_T_pt())[0] += w; */
// /* 	MolPt->get_W_T() += w / (*MolPt->get_R_T_pt())[0]; */
//
// /* 	// update the histogram */
// /* 	gsl_histogram_accumulate(hist,A,MolPt->get_W_T()); */
// /*       } */
// /*   } */
//
// /*   // rewieghting function that returns that adds the correct weight to the correct bin */
// /*   // see paper with Jianfeng, Eric and Ben. Always at lowest temperature */
// /*   void weighted_2d_histogram(const double& Ax,const double& Ay,  */
// /* 			     Molecule* MolPt, gsl_histogram2d* hist) */
// /*   { */
// /*     if(WithLearn) */
// /*       { */
// /* 	double V = *MolPt->potential_pt(); */
// /* 	std::vector<double> BetaWeight = *MolPt->get_beta_weight_pt(); */
//
// /* 	double W = 0.0; */
// /* 	unsigned nbeta = Beta.size(); */
//
// /* 	for(unsigned i=0;i<nbeta;i++) */
// /* 	  { */
// /* 	    W += GaussW[i] * BetaWeight[i] * exp(-Beta[i] * V); */
// /* 	  } */
//
// /* 	W *= BetaWeight[nbeta-1] * exp(-Beta[nbeta-1] * V); */
//
// /* 	// update the histogram */
// /* 	gsl_histogram2d_accumulate(hist,Ax,Ay,W); */
//
// /*       } */
// /*     else */
// /*       { */
// /* 	//senare */
// /*       } */
// /*   } */
// /* }; */
//
//
//
// /* /\****************************************************************************** */
// /*                            CT_BAOAB INTEGRATOR */
// /*  *****************************************************************************\/ */
// /* template <class SYSTEM> class CT_BAOAB : public Integrator<SYSTEM> */
// /* { */
// /*  private: */
// /*   double GammaBeta; */
// /*   double MassBeta; */
//
// /*   // B step, changes momentum - scales with beta variable! */
// /*   void B_step(Molecule* Mol_pt,const double& step_size) */
// /*   { */
// /*     // loop over the number of particles */
// /*     for(unsigned i=0;i<*Mol_pt->nparticle_pt();i++) */
// /*       { */
// /*   	// loop over the number of dimensions */
// /*   	for(unsigned j=0;j<*Mol_pt->dim_pt();j++) */
// /*   	  { */
// /*   	    *Mol_pt->particle_pt(i)->p_pt(j)=*Mol_pt->particle_pt(i)->p_pt(j) */
// /* 	      +(*Mol_pt->beta_pt()*(1/Mol_pt->kt()))*step_size*Integrator<SYSTEM>::Time_Step */
// /* 	      *(*Mol_pt->particle_pt(i)->f_pt(j)); */
// /*   	  } */
// /*       } */
// /*   }; */
//
// /*   // beta B step specific to this integrator */
// /*   void A_BETA_step(Molecule* Mol_pt, const double& stepsize) */
// /*   { */
// /*     *Mol_pt->beta_pt() = *Mol_pt->beta_pt() */
// /*       + stepsize * Integrator<SYSTEM>::Time_Step * (*Mol_pt->beta_mom_pt())/MassBeta; */
// /*   } */
//
// /*   // beta B step specific to this integrator */
// /*   void B_BETA_step(Molecule* Mol_pt, const double& stepsize) */
// /*   { */
// /*     // get partition function for current beta */
// /*     double weight = Integrator<SYSTEM>::System_pt-> */
// /*       domega_dbeta_over_omega(*Mol_pt->beta_pt(),*Mol_pt->dim_pt()); */
//
// /*     (*Mol_pt->beta_mom_pt()) = (*Mol_pt->beta_mom_pt())  */
// /*       - stepsize * Integrator<SYSTEM>::Time_Step * (Mol_pt->kt()) */
// /*       * ((*Mol_pt->potential_pt()) + weight); */
// /*   } */
//
// /*   // beta B step specific to this integrator */
// /*   void O_BETA_step(Molecule* Mol_pt, const double& stepsize) */
// /*   { */
// /*     // Calculate the c term */
// /*     double c=exp(-GammaBeta * stepsize * Integrator<SYSTEM>::Time_Step); */
//
// /*     (*Mol_pt->beta_mom_pt()) = c * (*Mol_pt->beta_mom_pt())  */
// /*       + sqrt((1-c*c) * (Mol_pt->kt())) */
// /*       * sqrt(MassBeta) * gsl_ran_gaussian(Integrator<SYSTEM>::r,1.0); */
// /*   } */
//
// /*  public: */
// /*   // Constructor */
// /*   CT_BAOAB(const double& time_step, */
// /* 	   const double& gamma, */
// /* 	   const double& gammabeta, */
// /* 	   const double& massbeta) */
// /*     { */
// /*       // Set the time step and Langevin Friction */
// /*       Integrator<SYSTEM>::Time_Step = time_step; */
// /*       Integrator<SYSTEM>::Gamma     = gamma; */
// /*       // set the beta values */
// /*       GammaBeta      = gammabeta; */
// /*       MassBeta       = massbeta; */
//
// /*       // create pointer to system */
// /*       Integrator<SYSTEM>::System_pt = new SYSTEM; */
//
// /*       // Generate a random seed based on system variables */
// /*       long unsigned seed; */
// /*       struct timeval tv; */
// /*       gettimeofday(&tv,0); */
// /*       seed = tv.tv_sec + tv.tv_usec; */
//
// /*       // Select and seed the random number generator  */
// /*       Integrator<SYSTEM>::r = gsl_rng_alloc (gsl_rng_mt19937); */
// /*       gsl_rng_set(Integrator<SYSTEM>::r,seed); */
// /*     }; */
//
// /*   // Integrator function */
// /*   void Integrate(Molecule* Mol_pt) */
// /*   { */
// // /*     // F is force from previous step */
// // /*     // The first steps at this force level */
// // /*     B_step(Mol_pt,0.5); */
// // /*     Integrator<SYSTEM>::A_step(Mol_pt,0.5); *x/
// // /*     B_BETA_step(Mol_pt,0.5); */
// // /*     A_BETA_step(Mol_pt,0.5); */
//
// // /*     O_BETA_step(Mol_pt,0.5); */
// // /*     Integrator<SYSTEM>::O_step(Mol_pt,1.0); */
//
// // /*     O_BETA_step(Mol_pt,0.5); */
// // /*     A_BETA_step(Mol_pt,0.5); */
//
// // /*     // potential evaluation */
// // /*     Integrator<SYSTEM>::System_pt->compute_potential(Mol_pt); */
//
// // /*     B_BETA_step(Mol_pt,0.5); */
// // /*     Integrator<SYSTEM>::A_step(Mol_pt,0.5); */
//
// // /*     // Force Solve */
// // /*     Integrator<SYSTEM>::System_pt->compute_force(Mol_pt); */
//
// // /*     // final B step for new force level */
// // /*     B_step(Mol_pt,0.5); */
// /*   } */
//
// /*   // Destructor */
// /*   ~CT_BAOAB()=delete; */
// /* }; */
//
// /* void simulated_tempering_switch(std::vector<Molecule*>& Molecules, */
// /* 				const std::vector<double>& beta, */
// /* 				const std::vector<double>& weight, */
// /* 				const std::vector<double>& Z_beta); */
//
//
// /****************************************************************************
//                             EULER MARUYAMA INTEGRATOR
// ****************************************************************************/
// template <class SYSTEM> class EulerMaruyama : public Integrator<SYSTEM>
// {
// private:
//
//   void step(Particle* PartPt, const double& kT)
//   {
//     // dereference
//     double* q = PartPt->q_pt(0);
//     double  f = *PartPt->f_pt(0);
//
//     double  R  = gsl_ran_gaussian(Integrator<SYSTEM>::r,1.0);
//     double h   = Integrator<SYSTEM>::Time_Step;
//
//     // senare loop over dimensions
//     *q += h * f + sqrt(2.0 * kT * h) * R;
//   }
//
// public:
//   EulerMaruyama(const double& time_step)
//   {
//     // initialise the integrator object - dummy friction gamma
//     Integrator<SYSTEM>::initialise_integrator(time_step,1.0,0.0);
//   };
//
//   void Integrate(Molecule* MolPt)
//   {
//     // determine the number of particles
//     unsigned nParticle = MolPt->nparticle();
//
//     // determine the temperature
//     double  kT = MolPt->kt();
//
//     // determine the dimension
//     double  Dim = MolPt->dim();
//
// #pragma omp single
//     Integrator<SYSTEM>::fleming_viot_count++;
//
//     // loop over all the particles
// #pragma omp for
//     for(unsigned i=0;i<nParticle;i++)
//       {
// 	// step forward in time
// 	step(MolPt->particle_pt(i),kT);
//
// 	// check boundary condition
// 	if(Integrator<SYSTEM>::WITH_FLEMING_VIOT
// 	   && (Integrator<SYSTEM>::fleming_viot_count % Integrator<SYSTEM>::fleming_viot_tau == 0))
// 	  {
// 	    int index = 0;
// 	    Integrator<SYSTEM>::fleming_viot_count=1;
//
// 	    while(index != -1)
// 	      {
// 		// get index
// 		index = Integrator<SYSTEM>::enforce_fleming_viot(MolPt,MolPt->particle_pt(i));
//
// 		// make a copy of the new particlex
// 		if(index != -1)
// 		  {
// 		    for(unsigned j=0;j<Dim;j++)
// 		      {
// 			*MolPt->particle_pt(i)->q_pt(j) = *MolPt->particle_pt(index)->q_pt(j);
// 			*MolPt->particle_pt(i)->f_pt(j) = *MolPt->particle_pt(index)->f_pt(j);
// 		      }
// 		  }
// 	      }
// 	  }
//
// 	// update the force
// 	Integrator<SYSTEM>::System_pt->compute_force(MolPt->particle_pt(i),Dim);
//       }
//   };
// };

#endif
