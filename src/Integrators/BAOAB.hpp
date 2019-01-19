#include "Integrators.hpp"
#include "System.hpp"

using namespace::std;

/******************************************************************************
                                 BAOAB Class
 *****************************************************************************/
class BAOAB : public Langevin
{
public:
    BAOAB(const double& beta, const double& gamma, const double& gamma_rot,
          const double& time_step, System* system_pt);
    // destructor
    ~BAOAB();
    // must have integrate
    void integrate(Molecule* molecule_pt);
    // set with npt grid
    void integrate_with_npt_grid(const double& a_x, const double& b_x,
                                 const double& b_y, const double& cut_off,
                                 Molecule* molecule_pt, const double& mass,
                                 const double& target_press,
                                 const double& gamma_rot);
    // integrate with simulated tempering
    void integrate_with_st(const double& tmin,
                           const double& tmax,
                           const double& n_temperatures,
                           const unsigned& mod_switch);
    // set which npt integration version
    void set_npt_integrator_version(const unsigned& version);

private:
    double Time_Step;
    unsigned Step;
    unsigned Npt_version;

    bool With_npt;
    bool With_st;

    void nvt_integration(Molecule* molecule_pt);
    void npt_integration(Molecule* molecule_pt);
};



// /******************************************************************************
//                                BAOAB INTEGRATOR
//
//  This file contains the -Langevin Integrator. This can be used with:
//  1. ISST
//  2. Rigid Body Dynamics
//  3. NPT thermostat
//
//  *****************************************************************************/
//
// #ifndef BAOAB_HPP
// #define BAOAB_HPP
//
// #include <vector>
// #include <omp.h>
//
// #include "Integrators.hpp"
// #include "Methods.hpp"
// #include "Thermostats.hpp"
//
//
// template <class SYSTEM> class BAOAB : public Integrator<SYSTEM>
// {
//  private:
//    // ------------------------- GRID VARIABLES ------------------------------ //
//    Grid* grid_pt;
//    bool with_grid_acceleration;
//
//    // ----------------------- NPT GRID VARIABLES ---------------------------- //
//    NptGrid* npt_grid_pt;
//    bool with_npt_grid_acceleration;
//    bool with_npt_method_nr;
//
//    NptThermostat* npt_thermostat;
//
//   // -------------------------- ISST VARIABLES ------------------------------ //
//   ISST* isst_pt;
//   bool with_isst_acceleration;
//
//   // ------------------ SIMULATED TEMPERING VARIABLES ----------------------- //
//   SimulatedTempering* simulated_tempering_pt;
//   bool with_simulated_tempering;
//
//   // ------------------------- PRIVATE FUNCTIONS ---------------------------- //
//   void pre_force_integration(Particle* PartPt, const unsigned& dim,
//                              const double& kT);
//
//    // integration strategy after force integration
//   void post_force_integration(Particle* PartPt, const unsigned& dim);
//
//   void A(Particle* PartPt, const unsigned& dim, const double& kT);
//
//   // local B step -- for experimental use
//   void B(Particle* PartPt, const unsigned& dim, const double& kT);
//
//   // local O step -- for experimental use
//   void O(Particle* PartPt, const unsigned& dim, const double& kT);
//
//  public:
//    //constructor
//    BAOAB(const double& time_step, const double& gamma, const double& gamma_rot,
//          const double& kT);
//
//    // destructor
//   ~BAOAB()=delete;
//
//   // integrate forward in time
//   void Integrate(Molecule* Mol_pt);
//   void standard_integration(Molecule* MolPt);
//   void npt_integration(Molecule* MolPt);
//    void npt_integration_v2(Molecule* MolPt);
//   void grid_integration(Molecule* MolPt);
//
//   // set isst method pointer
//   void set_ISST(Molecule* MolPt, const double& tmin, const double& tmax,
//                 const unsigned& nInt, const double& alpha);
//
//   // get the isst beta temperature
//   double get_isst_beta(const double& i);
//
//   // set all the beta weights
//   void set_isst_beta_weight(const std::vector<double>& weight,
//                             Molecule* molecule_pt);
//
//   // get all the weights for the observable
//   std::vector<double> get_isst_observable_weights(Molecule* molecule_pt);
//
//   // initialisor function for setting the grid
//   void set_grid(const double& a_x, const double& b_x, const double& b_y,
//                 const double& cut_off, Molecule* molecule_pt);
//
//   // initialisor function for setting the grid
//   void set_npt_grid(const double& a_x, const double& b_x, const double& b_y,
//                     const double& cut_off, Molecule* molecule_pt,
//                     const double& mass, const double& target_pressure,
//                     const double& npt_langevin_friction);
//
//   // initialisor function for setting the grid
//   void set_npt_initial(Molecule* molecule_pt,
//                        const char* initial_pos_filename,
//                        const char* initial_mom_filename,
//                        const char* initial_box_filename)
//
//   {
//       // read in the particles and the box
//       npt_grid_pt->add_file_initial_condition(molecule_pt, initial_pos_filename,
//                                               initial_mom_filename,
//                                               initial_box_filename);
//
//       // check the initial condition
//       npt_grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt,
//                                           molecule_pt);
//   }
//
//   // set simulated tempering
//   void set_simulated_tempering(const double& tmin,
//                                const double& tmax,
//                                const double& n_temperatures,
//                                const unsigned& mod_switch);
//
//   unsigned get_current_temperature_index();
//
//   double get_simulated_tempering_beta(const double& i);
//
//   void set_simulated_tempering_beta_weight(const double& weight,
//                                           const unsigned& i);
//
//   double get_simulated_tempering_temperature();
//
//   void update_pressure_temperature()
//   {
//       npt_grid_pt->update_pressure_temperature();
//   }
//
//   void set_npt_integration_method(const unsigned& nr)
//   {
//       if(nr == 1)
//         with_npt_method_nr = false;
//       else if(nr == 2)
//         with_npt_method_nr = true;
//       else
//       {
//         printf("\nERROR: uknown NPT integration!\n");
//         exit(-1);
//       }
//   }
//
//   void update_temperature()
//   {
//       grid_pt->update_temperature();
//   }
//
//   // returns the pressure in the npt grid
//   double get_npt_pressure();
//   double get_instant_npt_pressure()
//   {
//       return npt_grid_pt->get_instant_pressure();
//   }
//
//   // returns the temperature in the grid
//   double get_temperature();
//   double get_instant_temperature()
//   {
//       // holder for the temperature
//       double temp = 0.0;
//
//       // get the temperature from the correct pointer
//       if(with_grid_acceleration == true)
//       {
//           temp = grid_pt->get_instant_temperature();
//       }
//       else if(with_npt_grid_acceleration == true)
//       {
//           temp = npt_grid_pt->get_instant_temperature();
//       }
//
//       // return the temperature
//       return temp;
//   }
//
//   // returns the volumne of the grid
//   double get_npt_volume();
//   double get_instant_npt_volume()
//   {
//       return npt_grid_pt->get_instant_volume();
//   }
//
//   // returns vector with the box momentum
//   std::vector<double> get_npt_box_momentum()
//   {
//       std::vector<double> momentum(3, 0.0);
//
//       for(unsigned i=0; i<3; i++)
//             momentum[i] = *npt_grid_pt->Lp_pt(i);
//
//     return momentum;
//   }
//
//   // returns vector with the box momentum
//   std::vector<double> get_npt_box_coord()
//   {
//       std::vector<double> coordinate(3, 0.0);
//
//       for(unsigned i=0; i<3; i++)
//             coordinate[i] = *npt_grid_pt->L_pt(i);
//
//     return coordinate;
//   }
// };
//
// #endif
//
// // Constructor
// template<class SYSTEM> inline BAOAB<SYSTEM>::BAOAB(const double& time_step,
//                                                    const double& gamma,
//                                                    const double& gamma_rot,
//                                                    const double& kT)
// {
//   // initialise the integrator
//   Integrator<SYSTEM>::initialise_integrator(time_step, gamma, kT, 1.0,
//                                             gamma_rot);
//
//   // set default with ISST to false
//   with_isst_acceleration = false;
//
//   // set the default grid parameter
//   with_grid_acceleration = false;
//
//   // set the default npt grid parameter
//   with_npt_grid_acceleration = false;
//   with_npt_method_nr = false;
//
//   // set the default ST parameter
//   with_simulated_tempering = false;
// };
//
// /*****************************************************************************
//  //
//  //                           INTEGRATION FUNCTIONS
//  //
//  *****************************************************************************/
//
// // Integrator function
// template<class SYSTEM> inline void BAOAB<SYSTEM>::Integrate(Molecule* MolPt)
// {
//     // increase the timestep
// #pragma omp single
//     Integrator<SYSTEM>::TIMEINDEX++;
//
//     // determine how to integrate forward
//     if(with_grid_acceleration == true)
//         grid_integration(MolPt);
//     else if(with_npt_grid_acceleration == true)
//     {
//         if(with_npt_method_nr == false)
//             npt_integration(MolPt);
//         else if(with_npt_method_nr == true)
//             npt_integration_v2(MolPt);
//     }
//     else
//         standard_integration(MolPt);
//
//     // check different methods
//     if(with_simulated_tempering)
//     {
//         // update the temperature
//         simulated_tempering_pt->update_temperature(MolPt,
//                                                  Integrator<SYSTEM>::TIMEINDEX);
//
//         // update the o step parameters
//         double temperature = simulated_tempering_pt->get_temperature();
//         Integrator<SYSTEM>::update_o_step_parameters(temperature, 1.0);
//     }
// };
//
// // Integrator function
// template<class SYSTEM> inline void BAOAB<SYSTEM>::npt_integration(Molecule* molecule_pt)
// {
//   // determine the number of particles
//   unsigned nParticle = molecule_pt->nparticle();
//
//   // temperature and dimension
//   double kT = molecule_pt->kt();
//   unsigned rot_dim = 1;
//
//       //#pragma omp for schedule(static)
// #pragma omp barrier
// #pragma omp single
//     {
//         //-------------------------------------------------------------- B STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         #pragma omp simd
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           npt_thermostat->B(particle_k, npt_grid_pt, 0.5);
//
//           if(particle_k->rigid_body())
//             Integrator<SYSTEM>::BpPhi_step(particle_k, 0.5, rot_dim);
//         }
//
//         //------------------------------------------------------------- A1 STEP
//         npt_thermostat->A_one(npt_grid_pt, 0.5);
//
//         //------------------------------------------------------------- A2 STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         #pragma omp simd
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           npt_thermostat->A_two(particle_k, npt_grid_pt, 0.5);
//
//           if(particle_k->rigid_body())
//             Integrator<SYSTEM>::ApPhi_step(particle_k, 0.5, rot_dim);
//         }
//
//         //-------------------------------------------------------------- O STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         #pragma omp simd
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             npt_thermostat->O(particle_k, npt_grid_pt);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::OpPhi_step(particle_k, kT);
//         }
//
//         //------------------------------------------------------------- A2 STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         // All the reverse order A steps
//         #pragma omp simd
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::ApPhi_step(particle_k, 0.5, rot_dim);
//
//             npt_thermostat->A_two(particle_k, npt_grid_pt, 0.5);
//         }
//
//         //------------------------------------------------------------- A1 STEP
//         npt_thermostat->A_one(npt_grid_pt, 0.5);
//
//     } // end of omp single evaluation
//
//     // integrate the force forward (using parallelisation)
//     npt_grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt,
//                                         molecule_pt);
//
//     #pragma omp barrier
//     #pragma omp single
//     {
//         //-------------------------------------------------------------- B STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         // loop over particles
//         #pragma omp simd
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::BpPhi_step(particle_k, 0.5, rot_dim);
//
//             npt_thermostat->B(particle_k, npt_grid_pt, 0.5);
//         }
//     }
// };
//
// // Integrator function
// template<class SYSTEM> inline void BAOAB<SYSTEM>::npt_integration_v2(Molecule* molecule_pt)
// {
//   // determine the number of particles
//   unsigned nParticle = molecule_pt->nparticle();
//
//   // temperature and dimension
//   double kT = molecule_pt->kt();
//   unsigned rot_dim = 1;
//
//       //#pragma omp for schedule(static)
// #pragma omp barrier
// #pragma omp single
//     {
//         //-------------------------------------------------------------- B STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           npt_thermostat->B(particle_k, npt_grid_pt, 0.5);
//
//           if(particle_k->rigid_body())
//             Integrator<SYSTEM>::BpPhi_step(particle_k, 0.5, rot_dim);
//         }
//
//         //------------------------------------------------------------- A2 STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           npt_thermostat->A_two(particle_k, npt_grid_pt, 0.5);
//
//           if(particle_k->rigid_body())
//             Integrator<SYSTEM>::ApPhi_step(particle_k, 0.5, rot_dim);
//         }
//
//         //------------------------------------------------------------- A1 STEP
//         npt_thermostat->A_one(npt_grid_pt, 0.5);
//
//         //-------------------------------------------------------------- O STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             npt_thermostat->O(particle_k, npt_grid_pt);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::OpPhi_step(particle_k, kT);
//         }
//
//         //------------------------------------------------------------- A1 STEP
//         npt_thermostat->A_one(npt_grid_pt, 0.5);
//
//         //------------------------------------------------------------- A2 STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         // All the reverse order A steps
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::ApPhi_step(particle_k, 0.5, rot_dim);
//
//             npt_thermostat->A_two(particle_k, npt_grid_pt, 0.5);
//         }
//
//     } // end of omp single evaluation
//
//     // integrate the force forward (using parallelisation)
//     npt_grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt,
//                                         molecule_pt);
//
//     #pragma omp barrier
//     #pragma omp single
//     {
//         //-------------------------------------------------------------- B STEP
//         // set the box to integrate
//         npt_thermostat->set_box_to_integrate();
//
//         // loop over particles
//         for(unsigned k=0; k<nParticle; k++)
//         {
//             Particle* particle_k = molecule_pt->particle_pt(k);
//
//             if(particle_k->rigid_body())
//                 Integrator<SYSTEM>::BpPhi_step(particle_k, 0.5, rot_dim);
//
//             npt_thermostat->B(particle_k, npt_grid_pt, 0.5);
//         }
//     }
// };
//
//
// // Integrator function
// template<class SYSTEM> inline void BAOAB<SYSTEM>::grid_integration(Molecule* molecule_pt)
// {
//   // determine the number of particles
//   unsigned nParticle = molecule_pt->nparticle();
//   unsigned Dim = molecule_pt->dim();
//
//   // temperature and dimension
//   double kT = molecule_pt->kt();
//
//       //#pragma omp for schedule(static)
// #pragma omp barrier
// #pragma omp single
//     {
//         // All the pre force steps
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           // perform all the post force integration
//           pre_force_integration(particle_k, Dim, kT);
//         }
//     }
//
//     // integrate the force forward (using parallelisation)
//     grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt,
//                                     molecule_pt);
//
//     #pragma omp barrier
//     #pragma omp single
//     {
//         // All the post force steps
//         for(unsigned k=0; k<nParticle; k++)
//         {
//           Particle* particle_k = molecule_pt->particle_pt(k);
//
//           // perform all the post force integration
//           post_force_integration(particle_k, Dim);
//         }
//     }
// };
//
// // Integrator function
// template<class SYSTEM> inline void BAOAB<SYSTEM>::standard_integration(Molecule* MolPt)
// {
//   // determine the number of particles
//   unsigned nParticle = MolPt->nparticle();
//
//   // temperature and dimension
//   double kT = MolPt->kt();
//   unsigned Dim = MolPt->dim();
//
//   // preforce integration
//   #pragma omp barrier
//   #pragma omp single
//   for(unsigned i=0;i<nParticle;i++)
//   {
//     // F force is from the previous step
//     Particle* PartPt = MolPt->particle_pt(i);
//
//     // perform all the pre force steps
//     pre_force_integration(PartPt, Dim, kT);
//   }
//
//   // Force Solve
//   if(with_grid_acceleration)
//   {
//       // this is not safe for parallelisation just yet
//       #pragma omp barrier
//       grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt, MolPt);
//   }
//   else
//       Integrator<SYSTEM>::System_pt->compute_force(MolPt);
//
//   // scale the force according ISST
//   #pragma omp barrier
//   #pragma omp single
//   if(with_isst_acceleration)
//     {
//         // compute the potential -- senare make sure that this is allways called
//         // in the force calculation
//     	Integrator<SYSTEM>::System_pt->compute_potential(MolPt);
//
//         // apply the force rescaling
//         isst_pt->apply_force_rescaling(MolPt);
//     }
//
//   // post force integrations
//   #pragma omp barrier
//   #pragma omp single
//   for(unsigned i=0;i<nParticle;i++)
//     {
//       // F force is from the previous step
//       Particle* PartPt = MolPt->particle_pt(i);
//
//       // perform all the post force integration
//       post_force_integration(PartPt, Dim);
//     }
// };
//
//
// // Integrator function particles pre_force
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::pre_force_integration(Particle* PartPt,
//                                           const unsigned& dim,
//                                           const double& kT)
// {
//   // step forward BAOA
//   Integrator<SYSTEM>::B_step(PartPt, 0.5, dim);
//   Integrator<SYSTEM>::A_step(PartPt, 0.5, dim);
//   Integrator<SYSTEM>::O_step(PartPt, 1.0, kT, dim);
//   Integrator<SYSTEM>::A_step(PartPt, 0.5, dim);
//
//   // if particle also rotates, integrate forward
//   // the rotation too.
//   if(PartPt->rigid_body())
//     {
//       Integrator<SYSTEM>::BpPhi_step(PartPt, 0.5, 1);
//       Integrator<SYSTEM>::ApPhi_step(PartPt, 0.5, 1);
//       Integrator<SYSTEM>::OpPhi_step(PartPt, kT);
//       Integrator<SYSTEM>::ApPhi_step(PartPt, 0.5, 1);
//     }
// };
//
// // post force integration
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::post_force_integration(Particle* PartPt,
//                                            const unsigned& dim)
// {
//     // step forward B
//     Integrator<SYSTEM>::B_step(PartPt, 0.5, dim);
//
//     // if particles rotates
//     if(PartPt->rigid_body())
//     {
//         Integrator<SYSTEM>::BpPhi_step(PartPt, 0.5, 1);
//     }
// };
//
//
// // A step
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::A(Particle* PartPt,
// 		      const unsigned& dim,
// 		      const double& kT)
// {
//   Integrator<SYSTEM>::A_step(PartPt,0.5,dim);
//
//   if(PartPt->rigid_body() != false)
//     Integrator<SYSTEM>::ApPhi_step(PartPt,0.5);
// };
//
// // B step
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::B(Particle* PartPt,
// 		      const unsigned& dim,
// 		      const double& kT)
// {
//   Integrator<SYSTEM>::B_step(PartPt,0.5,dim);
//
//   if(PartPt->rigid_body() != false)
//     Integrator<SYSTEM>::BpPhi_step(PartPt,0.5);
// };
//
// // O step
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::O(Particle* PartPt,
// 		      const unsigned& dim,
// 		      const double& kT)
// {
//   Integrator<SYSTEM>::O_step(PartPt,
// 			     1.0,
// 			     kT,
// 			     dim);
//
//   if(PartPt->rigid_body() != false)
//     Integrator<SYSTEM>::OpPhi_step(PartPt,kT);
//
// };
//
// /*****************************************************************************
//  //
//  //                               ISST FUNCTIONS
//  //
//  *****************************************************************************/
// // initialise ISST
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::set_ISST(Molecule* MolPt,
//                              const double& tmin,
// 			                 const double& tmax,
// 			                 const unsigned& nInt,
// 			                 const double& tau)
// {
//   // scale the timestep for of the learning rate
//   double time_step = Integrator<SYSTEM>::Time_Step;
//
//   // initialise the class of methods
//   isst_pt = new ISST(MolPt, tmin, tmax, nInt, time_step, tau);
//
//   // set boolean for With_ISST to true
//   with_isst_acceleration = true;
//
//   // get the method info
//   isst_pt->get_method_info();
// };
//
// template<class SYSTEM> inline
// double BAOAB<SYSTEM>::get_isst_beta(const double& i)
// {
//     return isst_pt->get_beta(i);
// }
//
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::set_isst_beta_weight(const std::vector<double>& weight,
//                                          Molecule* molecule_pt)
// {
//     isst_pt->set_beta_weight(weight, molecule_pt);
// }
//
// template<class SYSTEM> inline
// std::vector<double> BAOAB<SYSTEM>::get_isst_observable_weights(Molecule* molecule_pt)
// {
//     return isst_pt->get_observable_weights(molecule_pt);
// }
//
// /*****************************************************************************
//  //
//  //                           GRID ACCELERATION
//  //
//  *****************************************************************************/
//  // initialise ISST
//  template<class SYSTEM> inline
//  void BAOAB<SYSTEM>::set_grid(const double& a_x, const double& b_x,
//                               const double& b_y, const double& cut_off,
//                               Molecule* molecule_pt)
// {
//   // build a new grid
//   grid_pt = new Grid(a_x, b_x, b_y, cut_off, molecule_pt);
//
//   // set boolean with grid to true
//   with_grid_acceleration = true;
//
//   // tell user we are running initial condition test
//   printf("Integrator is testing initial condition:");
//
//   // check the initial condition
//   grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt, molecule_pt);
//
//   // let user know that the initial condiont was accepted
//   printf("...accepted!\n");
// };
//
// /*****************************************************************************
//  //
//  //                         NPT GRID ACCELERATION
//  //
//  *****************************************************************************/
//  // initialise ISST
//  template<class SYSTEM> inline
//  void BAOAB<SYSTEM>::set_npt_grid(const double& a_x, const double& b_x, const double& b_y,
//                    const double& cut_off, Molecule* molecule_pt,
//                    const double& mass, const double& target_pressure,
//                    const double& npt_langevin_friction)
// {
//   // build a new grid
//   npt_grid_pt = new NptGrid(a_x, b_x, b_y, cut_off, molecule_pt, mass,
//                             target_pressure);
//
//   // set boolean with grid to true
//   with_npt_grid_acceleration = true;
//
//   // tell user we are running initial condition test
//   printf("Integrator is testing initial condition:");
//
//   // check the initial condition
//   npt_grid_pt->update_particle_forces(Integrator<SYSTEM>::System_pt,
//                                       molecule_pt);
//
//   // update the box forces
//   // npt_grid_pt->update_box_force();
//
//   // let user know that the initial condiont was accepted
//   printf("...accepted!\n");
//
//   // build new thermostat
//   double kt = molecule_pt->kt();
//   double time_step = Integrator<SYSTEM>::Time_Step;
//
//   npt_thermostat = new NptThermostat(target_pressure, kt,
//                                      Integrator<SYSTEM>::Gamma, 1.0,
//                                      npt_langevin_friction, time_step);
// };
//
// template<class SYSTEM> inline
//  double BAOAB<SYSTEM>::get_npt_pressure()
//  {
//    return npt_grid_pt->get_pressure();
//  }
//
//  template<class SYSTEM> inline
//   double BAOAB<SYSTEM>::get_temperature()
//   {
//       // holder for the temperature
//       double temp = 0.0;
//
//       // get the temperature from the correct pointer
//       if(with_grid_acceleration == true)
//         temp = grid_pt->get_temperature();
//       else if(with_npt_grid_acceleration == true)
//         temp = npt_grid_pt->get_temperature();
//
//       // return the temperature
//       return temp;
//   }
//
//   template<class SYSTEM> inline
//    double BAOAB<SYSTEM>::get_npt_volume()
//    {
//      return npt_grid_pt->get_volume();
//    }
//
// /*****************************************************************************
// //
// //                           SIMULATED TEMPERING
// //
// *****************************************************************************/
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::set_simulated_tempering(const double& tmin,
//                                             const double& tmax,
//                                             const double& n_temperatures,
//                                             const unsigned& mod_switch)
// {
//     // set boolean of simulated tempering
//     with_simulated_tempering = true;
//
//     // initialise the simulated tempering class
//     simulated_tempering_pt = new SimulatedTempering(tmin, tmax, n_temperatures,
//                                                     mod_switch);
//
//     // get the method info
//     simulated_tempering_pt->get_method_info();
// };
//
// template<class SYSTEM> inline
// unsigned BAOAB<SYSTEM>::get_current_temperature_index()
// {
//     return simulated_tempering_pt->get_temperature_index();
// }
//
// template<class SYSTEM> inline
// double BAOAB<SYSTEM>::get_simulated_tempering_beta(const double& i)
// {
//     return simulated_tempering_pt->get_beta(i);
// }
//
// template<class SYSTEM> inline
// double BAOAB<SYSTEM>::get_simulated_tempering_temperature()
// {
//     return simulated_tempering_pt->get_temperature();
// }
//
// template<class SYSTEM> inline
// void BAOAB<SYSTEM>::set_simulated_tempering_beta_weight(const double& weight,
//                                                        const unsigned& i)
// {
//     simulated_tempering_pt->set_beta_weight(weight, i);
// }
