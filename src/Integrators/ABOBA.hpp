// /******************************************************************************
//                                ABOBA INTEGRATOR
//
//  This file contains the ABOBA-Langevin Integrator. This can be used with:
//  1. ISST
//  2. Rigid Body Dynamics
//  3. NPT thermostat
//
// *****************************************************************************/
//
// #ifndef ABOBA_HPP
// #define ABOBA_HPP
//
// #include <vector>
// #include <omp.h>
//
// #include "Integrators.hpp"
// #include "Methods.hpp"
// #include "Thermostats.hpp"
//
//
// template <class SYSTEM> class ABOBA : public Integrator<SYSTEM>
// {
// private:
//
//  //  std::vector<double> Beta;                    // inverse temps
//  //  std::vector<double> GaussW;                  // Gauss Weights
//  //
//  //  // const gsl_rng_type* T;                       // ISST Random Seed
//  //  gsl_rng*            r_ISST;                  // ISST Random Seed
//  //
//  //  NPT* THERMOSTAT_NPT;                  // Constant pressure and temperature
//  //  bool     With_NPT;                           // Logical NPT
//  //
//  //  Methods* METHOD_ISST;                        // ISST Method pointer
//  //  bool     With_ISST;                          // Logical ISST
//  //
//  //  Grid* GRID_PT;    // pointer for accelerating the calucluation in parallel
//  //
//  // public:
//
//   ABOBA
//   (const double& time_step,
//    const double& gamma,
//    const double& kT);
//
//   ~ABOBA()=delete;                             // Destructor
//
//   void Integrate                               // Integration: Molecule
//   (Molecule* Mol_pt);
//
//   // void Integrate_NPT
//   // (Molecule* MolPt);
//   //
//   // void pre_force_integration                   // Integration: Particle
//   // (Particle*       PartPt,                     // pre force
//   //  const unsigned& dim,
//   //  const double&   kT);
//   //
//   // void post_force_integration                  // Integration: Particle
//   // (Particle*       PartPt,                     // post force
//   //  const unsigned& dim);
//   //
//   // void A                                       // Integration: Particle
//   // (Particle*       PartPt,                     // A step
//   //  const unsigned& dim,
//   //  const double&   kT);
//   //
//   // void B                                       // Integration: Particle
//   // (Particle*       PartPt,                     // B step
//   //  const unsigned& dim,
//   //  const double&   kT);
//   //
//   // void O                                       // Integration: Particle
//   // (Particle*       PartPt,                     // O step
//   //  const unsigned& dim,
//   //  const double&   kT);
//   //
//   // void set_ISST                                // ISST initialiser
//   // (Molecule*       MolPt,
//   //  const double&   tmin,
//   //  const double&   tmax,
//   //  const unsigned& nInt,
//   //  const double&   alpha);
//   //
//   // ISST* ISST_method_pt();
//   //
//   // void set_NPT                                 // NPT initialiser
//   // (Molecule* MolPt,
//   //  const double& L,
//   //  const double& pressure,
//   //  const double& CutOff,
//   //  const double& mu_p,
//   //  const double& eta);
// };
//
// #endif
//
// // Constructor
// template<class SYSTEM> inline
// ABOBA<SYSTEM>::ABOBA(const double& time_step,const double& gamma,
// 		     const double& kT)
// {
//   // initialise the integrator
//   Integrator<SYSTEM>::initialise_integrator(time_step,gamma,kT,1.0);
//
//   // // set default with ISST to false
//   // With_ISST = false;
//   //
//   // // set default with NPT thermostat to false
//   // With_NPT = false;
// };
//
// /*****************************************************************************
//  //
//  //                           INTEGRATION FUNCTIONS
//  //
//  *****************************************************************************/
// //
// // // Integrator function
// // template<class SYSTEM> inline
// // void ABOBA<SYSTEM>::Integrate(Molecule* MolPt)
// // {
// //   // increase the timestep
// // #pragma omp single
// //   Integrator<SYSTEM>::TIMEINDEX++;
// //
// //   // determine the number of particles
// //   unsigned nParticle = MolPt->nparticle();
// //
// //   // temperature and dimension
// //   double kT    = MolPt->kt();
// //   unsigned Dim = MolPt->dim();
// //
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nParticle;i++)
// //     {
// //       // F force is from the previous step
// //       Particle* PartPt = MolPt->particle_pt(i);
// //
// //       // compute the first step
// //       Integrator<SYSTEM>::A_step(PartPt,0.5,Dim);
// //     }
// //
// //   // Force Solve
// //   Integrator<SYSTEM>::System_pt->compute_force(MolPt);
// //
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nParticle;i++)
// //     {
// //       // F force is from the previous step
// //       Particle* PartPt = MolPt->particle_pt(i);
// //
// //       // last four steps
// //       Integrator<SYSTEM>::B_step(PartPt,0.5,Dim);
// //       Integrator<SYSTEM>::O_step(PartPt,1.0,kT,Dim);
// //       Integrator<SYSTEM>::B_step(PartPt,0.5,Dim);
// //       Integrator<SYSTEM>::A_step(PartPt,0.5,Dim);
// //     }
// //
// //   // check Fleming viot conditions
// //   if(Integrator<SYSTEM>::WITH_FLEMING_VIOT)
// //     {
// //       // enforce the constraint
// //       Integrator<SYSTEM>::enforce_fleming_viot(MolPt);
// //     }
// //
// //   // check if we run with collecting trans_stats
// //   if(Integrator<SYSTEM>::WITH_TRANS_STATS)
// //     {
// //       Integrator<SYSTEM>::collect_trans_stat(MolPt);
// //     }
// // };
// //
// // // // Integrator function particles pre_force
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::pre_force_integration(Particle* PartPt,
// // // 					  const unsigned& dim,
// // // 					  const double& kT)
// // // {
// // //   printf("IMPLEMETN THIS\n");
// // //   exit(-1);
// // //   // Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
// // //   // Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
// // //   // Integrator<SYSTEM>::O_step(PartPt,1.0,dim,kT);
// // //   // Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
// //
// // //   // if(PartPt->rigid_body() != false)
// // //   //   {
// // //   //     // Integrate angular momentum forward
// // //   //     Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
// //
// // //   //     // Integrate rotation forward
// // //   //     Integrator<SYSTEM>::ApPhi_step(PartPt,1.0);
// // //   //   }
// //
// // // };
// //
// // // // post force integration
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::post_force_integration(Particle* PartPt,
// // // 					    const unsigned& dim)
// // // {
// // //   printf("IMPLEMETN THIS\n");
// // //   exit(-1);
// //
// //
// // //   // Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
// //
// // //   // if(PartPt->rigid_body() != false)
// // //   //   {
// // //   //     // integrate angular momentum forward
// // //   //     Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
// // //   //   }
// // // };
// //
// //
// // // // A step
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::A(Particle* PartPt,
// // // 		      const unsigned& dim,
// // // 		      const double& kT)
// // // {
// // //   Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
// //
// // //   if(PartPt->rigid_body() != false)
// // //     Integrator<SYSTEM>::ApPhi_step(PartPt,0.5);
// // // };
// //
// // // // B step
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::B(Particle* PartPt,
// // // 		      const unsigned& dim,
// // // 		      const double& kT)
// // // {
// // //   Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
// //
// // //   if(PartPt->rigid_body() != false)
// // //     Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
// // // };
// //
// // // // O step
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::O(Particle* PartPt,
// // // 		      const unsigned& dim,
// // // 		      const double& kT)
// // // {
// // //   Integrator<SYSTEM>::O_step(PartPt,1.0,dim,kT);
// //
// // //   if(PartPt->rigid_body() != false)
// // //     Integrator<SYSTEM>::OpPhi_step(PartPt,0.5,kT);
// //
// // // };
// //
// //
// // // /*****************************************************************************
// // //  //
// // //  //                               ISST FUNCTIONS
// // //  //
// // //  *****************************************************************************/
// //
// // // // return the pointer of the isst methdod
// // // template<class SYSTEM> inline
// // // ISST* ABOBA<SYSTEM>::ISST_method_pt()
// // // {
// // //   return dynamic_cast<ISST*>(METHOD_ISST);
// // // }
// //
// // // // initialise ISST
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::set_ISST(Molecule* MolPt,
// // // 			     const double& tmin,
// // // 			     const double& tmax,
// // // 			     const unsigned& nInt,
// // // 			     const double& alpha)
// // // {
// // //   // scale the timestep for of the learning rate
// // //   double time_step = alpha * Integrator<SYSTEM>::Time_Step;
// //
// // //   // initialise the class of methods
// // //   METHOD_ISST = new ISST(MolPt,tmin,tmax,nInt,time_step,alpha);
// //
// // //   // set boolean for With_ISST to true
// // //   With_ISST = true;
// // // };
// //
// // // /*****************************************************************************
// // //  //
// // //  //                           NPT FUNCTIONS
// // //  //
// // //  *****************************************************************************/
// //
// // // // void set thermostat
// // // template<class SYSTEM> inline
// // // void ABOBA<SYSTEM>::set_NPT(Molecule* MolPt,
// // 			    const double& L,
// // 			    const double& pressure,
// // 			    const double& CutOff,
// // 			    const double& mu_p,
// // 			    const double& eta)
// // {
// //   // set boolean flag
// //   With_NPT = true;
//
// //   // set the Cell method
// //   GRID_PT = new Grid(L,L,CutOff);
//
// //   // set the NPT thermostat in square box
// //   THERMOSTAT_NPT = new NPT(L,L,pressure,mu_p,eta,Integrator<SYSTEM>::Time_Step);
//
// //   // update the volume of the grid
// //   // this is so that we set all the particles in the grid
// //   THERMOSTAT_NPT->update_volume(MolPt,GRID_PT);
//
// //   //update the particles on the grid
// //   GRID_PT->update_particles_on_grid(MolPt);
//
// //   // get the number of boxes on the grid
// //   double nCells = GRID_PT->get_tot_nCells();
//
// //   // update the particles forces
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       GRID_PT->update_particle_forces(Integrator<SYSTEM>::System_pt,MolPt,i);
// //     }
// // };
//
// // // rescale step
// // // template<class SYSTEM> inline
// // // void BAOAB<SYSTEM>::rescale_q(Particle* part,
// // // 			      const double& V_prime,
// // // 			      const double& V,
// // // 			      const unsigned& dim)
// // // {
// // //   for(unsigned i=0;i<dim;i++)
// // //     {
// // //       *part->q_pt(i) *= pow(V_prime/V, 1.0/(double)dim);
// // //     }
// // // };
//
//
// // // // A prime step
// // // template<class SYSTEM> inline
// // // void BAOAB<SYSTEM>::A_prime(Molecule* MolPt,
// // // 			    const double& step_size)
// // // {
// // //   // dereference the old volume for rescaling
// // //   double V   = Volume;
// // //   unsigned N = MolPt->nparticle();
//
// // //   unsigned dim = MolPt->dim();
//
// // //   // inegrate Volume
// // //   Volume += step_size * Integrator<SYSTEM>::Time_Step * Xi/mu_P;
//
// // //   // rescale the positions
// // //   for(unsigned i=0;i<N;i++)
// // //     {
// // //       rescale_q(MolPt->particle_pt(i),Volume,V,dim);
// // //     }
// // // };
//
// // // // B step - particle
// // // template<class SYSTEM> inline
// // // void BAOAB<SYSTEM>::B_prime(Molecule* MolPt,
// // // 			    const double& step_size)
// // // {
// // //   // parameters
// // //   double q = 0.0;
// // //   double f = 0.0;
//
// // //   double N     = MolPt->nparticle();
// // //   unsigned dim = MolPt->dim();
//
// // //   double p = 0.0;
// // //   double m = 0.0;
//
// // //   double q_f = 0.0;
// // //   double p_m_p = 0.0;
//
// // //   for(unsigned i=0;i<N;i++)
// // //     {
// // //       for(unsigned j=0;j<dim;j++)
// // // 	{
// // // 	  q = *MolPt->particle_pt(i)->q_pt(j);
// // // 	  f = *MolPt->particle_pt(i)->f_pt(j);
//
// // // 	  p = *MolPt->particle_pt(i)->p_pt(j);
// // // 	  m = *MolPt->particle_pt(i)->m_pt(j);
//
// // // 	  q_f   += q * f;
// // // 	  p_m_p += p * p / m;
// // // 	}
// // //     }
//
// // //   // calculate new xi
// // //   Xi += step_size * Integrator<SYSTEM>::Time_Step * (1.0/(Volume * (double)dim) * (q_f + p_m_p) - Pressure);
// // // };
//
// // // // O step - particle
// // // template<class SYSTEM> inline
// // // void BAOAB<SYSTEM>::O_prime(const double& kT,
// // // 			    const double& step_size)
// // // {
// // //   double c   = exp(-Eta * step_size * Integrator<SYSTEM>::Time_Step);
//
// // //   Xi = c * Xi + sqrt((1 - c*c) * kT)
// // //     * sqrt(mu_P) *gsl_ran_gaussian(r_xi,1.0);
// // // };
//
// // // // return the box length
// // // template<class SYSTEM> inline
// // // double BAOAB<SYSTEM>::L(const unsigned& dim)
// // // {
// // //   return pow(Volume,1.0/(double)dim);
// // // };
//
// // template<class SYSTEM> inline
// // void ABOBA<SYSTEM>::Integrate_NPT(Molecule* MolPt)
// // {
// //   // senare below is the BAOAB integrator
// //   printf("IMPLEMETN THIS\n");
// //   exit(-1);
//
//
// //   // dereference the number of cells
// //   unsigned nCells = GRID_PT->get_tot_nCells();
//
// //   // dereference dimension and temperature
// //   double Temp = MolPt->kt();
// //   double Dim  = MolPt->dim();
//
// //   //----------------------------------- B STEP --------------------------------------//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       // number of particles in cell
// //       unsigned nPartCell = GRID_PT->get_box(i)->get_nparticles();
//
// //       // particle holder
// //       Particle* Part = NULL;
//
// //       for(unsigned k=0;k<nPartCell;k++)
// // 	{
// // 	  // get the particle and integrate forward using B
// // 	  Part = GRID_PT->get_box(i)->get_particle_pt(k);
// // 	  this->B(Part,Dim,Temp);
//
// // 	  // check periodic boundary condition
// // 	  GRID_PT->check_periodic_BC(Part);
// // 	}
// //     }
//
// //   //senare this is slowing things down!
// // #pragma omp barrier
// // #pragma omp single
// //   GRID_PT->update_particles_on_grid(MolPt);
//
// //   //----------------------------------- B' STEP --------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   THERMOSTAT_NPT->B(MolPt,0.5);
//
// //   //----------------------------------- A STEP --------------------------------------//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       // number of particles on
// //       unsigned nPartCell = GRID_PT->get_box(i)->get_nparticles();
//
// //       // particle holder
// //       Particle* Part = NULL;
//
// //       for(unsigned k=0;k<nPartCell;k++)
// // 	{
// // 	  Part = GRID_PT->get_box(i)->get_particle_pt(k);
// // 	  this->A(Part,Dim,Temp);
// // 	}
// //     }
//
// //   //----------------------------------- A' STEP --------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   {
// //     // integrate forward in time
// //     THERMOSTAT_NPT->A(MolPt,0.5);
//
// //     // update the Grid with the current volume
// //     THERMOSTAT_NPT->update_volume(MolPt,GRID_PT);
// //   }
//
// //   // update all the threads nCell values
// //   nCells = GRID_PT->get_tot_nCells();
//
// //   //----------------------------------- O STEP --------------------------------------//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       // number of particles on
// //       unsigned nPartCell= GRID_PT->get_box(i)->get_nparticles();
//
// //       // particle holder
// //       Particle* Part = NULL;
//
// //       for(unsigned k=0;k<nPartCell;k++)
// // 	{
// // 	  Part = GRID_PT->get_box(i)->get_particle_pt(k);
// // 	  this->O(Part,Dim,Temp);
// // 	}
// //     }
//
// //   //----------------------------------- O' STEP --------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   THERMOSTAT_NPT->O(Temp,1.0);
//
// //   //----------------------------------- A' STEP --------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   {
// //     // integrate forward in time
// //     THERMOSTAT_NPT->A(MolPt,0.5);
//
// //     // update the Grid with the current volume
// //     THERMOSTAT_NPT->update_volume(MolPt,GRID_PT);
// //   }
//
// //   // update all the threads nCell values
// //   nCells = GRID_PT->get_tot_nCells();
//
// //   //----------------------------------- A STEP --------------------------------------//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       // number of particles on
// //       unsigned nPartCell = GRID_PT->get_box(i)->get_nparticles();
//
// //       // particle holder
// //       Particle* Part = NULL;
//
// //       for(unsigned k=0;k<nPartCell;k++)
// // 	{
// // 	  Part = GRID_PT->get_box(i)->get_particle_pt(k);
// // 	  this->A(Part,Dim,Temp);
// // 	}
// //     }
//
// //   //----------------------------- FORCE EVALUATION ------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   GRID_PT->clear_particle_forces(MolPt);
//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     GRID_PT->update_particle_forces(Integrator<SYSTEM>::System_pt,MolPt,i);
//
// //   //----------------------------------- B' STEP --------------------------------------//
// // #pragma omp barrier
// // #pragma omp single
// //   THERMOSTAT_NPT->B(MolPt,0.5);
//
// //   //----------------------------------- B STEP --------------------------------------//
// // #pragma omp for schedule(static)
// //   for(unsigned i=0;i<nCells;i++)
// //     {
// //       // number of particles in cell
// //       unsigned nPartCell = GRID_PT->get_box(i)->get_nparticles();
//
// //       // particle holder
// //       Particle* Part = NULL;
//
// //       for(unsigned k=0;k<nPartCell;k++)
// // 	{
// // 	  // get the particle and integrate forward using B
// // 	  Part = GRID_PT->get_box(i)->get_particle_pt(k);
// // 	  this->B(Part,Dim,Temp);
//
// // 	  // check periodic boundary condition
// // 	  GRID_PT->check_periodic_BC(Part);
// // 	}
// //     }
//
// //   //senare this is slowing things down!
// // #pragma omp barrier
// // #pragma omp single
// //   GRID_PT->update_particles_on_grid(MolPt);
// // };
