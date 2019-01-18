// /******************************************************************************
//                                   Thermostats
//
//  This file contains.
//  1. Generic Thermostat class
//  2. 2D NPT thermostat
//
//  *****************************************************************************/
//
// #ifndef THERMOSTAT_HPP
// #define THERMOSTAT_HPP
//
// #include <vector>
// #include <omp.h>
// #include <stdio.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <sys/time.h>
// #include <math.h>
//
// #include "Molecules.hpp"
// #include "Cells.hpp"
//
// class Thermostat
// {
// public:
//   // empty constructor
//   Thermostat(){};
//
//   // empty destructor
//   ~Thermostat(){};
// };
//
// // The NPT thermostat, which derives from the general
// // thermostat class can be used with the Integrator class
// // to enforce:
// // N - Constant number of particles
// // P - Constant pressure
// // T - Constant Temperature
//
// class NptThermostat : public Thermostat
// {
// private:
//   // holder for the target pressure
//   double target_pressure;
//
//   // holder for the number of coordinates to
//   // integrate forward
//   unsigned number_of_grid_coordinates;
//
//   // holder for the dimension of the particles
//   unsigned dimension;
//
//   // holder for the time_step
//   double time_step;
//
//   // particle o step holders
//   double o_part_step_c;
//   double o_part_step_zeta;
//
//   // box o step holders
//   double o_box_step_c;
//   double o_box_step_zeta;
//
//   // holders for three different random number
//   std::vector<gsl_rng*> box_generator;
//   std::vector<gsl_rng*> part_generator;
//
//   // boolean for doing box integration step
//   bool box_integration;
//
// public:
//   // constructor
//   NptThermostat(const double& pressure, const double& temp,
//                 const double& part_langevin_friction,
//                 const double& o_step_size,
//                 const double& box_langevin_friction,
//                 const double& integrator_time_step)
//   {
//       // set the dimension of the box
//       dimension = 2;
//
//       // set the target pressure and volume
//       target_pressure = pressure;
//
//       // set the timestep of the thermostat
//       time_step = integrator_time_step;
//
//       // number of box coordinates
//       number_of_grid_coordinates = 3;
//
//       // set the box to integrate
//       box_integration = false;
//
//       // a seed for the random number generators
//       long unsigned seed;
//       struct timeval tv;
//       gettimeofday(&tv,0);
//       seed = tv.tv_sec + tv.tv_usec + 19.0 * omp_get_thread_num();
//
//       // set the particle generators
//       for(unsigned i=0; i<dimension; i++)
//       {
//           part_generator.push_back(gsl_rng_alloc(gsl_rng_mt19937));
//           gsl_rng_set(part_generator[i], seed + i * 3);
//       }
//
//       // set the three box generators
//       for(unsigned i=0; i<number_of_grid_coordinates; i++)
//       {
//           box_generator.push_back(gsl_rng_alloc(gsl_rng_mt19937));
//           gsl_rng_set(box_generator[i], seed + i * 4);
//       }
//
//       // set up the particle o-step constatns
//       o_part_step_c = exp(-part_langevin_friction * o_step_size * time_step);
//       o_part_step_zeta = sqrt((1.0 - o_part_step_c * o_part_step_c) * temp);
//
//       // set up the box o-step constatns
//       o_box_step_c = exp(-box_langevin_friction * o_step_size * time_step);
//       o_box_step_zeta = sqrt((1.0 - o_box_step_c * o_box_step_c) * temp);
//   }
//
//   // destructor
//   ~NptThermostat()
//   {
//     // free the random number generators
//     for(unsigned i=0; i<dimension; i++)
//       gsl_rng_free(part_generator[i]);
//     part_generator.clear();
//
//     for(unsigned i=0; i<number_of_grid_coordinates; i++)
//       gsl_rng_free(box_generator[i]);
//     box_generator.clear();
//   }
//
//   // Particle A step
//   void A_two(Particle* particle_pt, NptGrid* npt_grid_pt, const double& step_size);
//
//   // Particle B step
//   void B(Particle* particle_pt, NptGrid* npt_grid_pt, const double& step_size);
//
//   // Particle O step
//   void O(Particle* particle_pt, NptGrid* npt_grid_pt);
//
//   // set box to integrate forward in time
//   void set_box_to_integrate(){box_integration = true;}
//
//   // Grid A step
//   void A_one(NptGrid* npt_grid_pt, const double& step_size)
//   {
//     // make references
//     double* L = NULL;
//     double Lp = 0.0;
//     double mass = npt_grid_pt->get_mass();
//
//     // old position needed for rescaling
//     std::vector<double> L_old(3, 0.0);
//
//     for(unsigned i=0; i<number_of_grid_coordinates; i++)
//     {
//       // derefernce the grid coordinates
//       L = npt_grid_pt->L_pt(i);
//       Lp = *npt_grid_pt->Lp_pt(i);
//
//       // save the old position
//       L_old[i] = *L;
//
//       // update the coorinate
//       *L += step_size * time_step * Lp / mass;
//     }
//
//     // need to rescale the particles such that they are
//     // in the same position with repsect to the cell
//     npt_grid_pt->enforce_constant_relative_particle_pos(L_old);
//   }
//   //
//   // void B(NptGrid* npt_grid_pt, const double& step_size)
//   // {
//   //   // helpers
//   //   double* Lp = NULL;
//   //   double f = 0.0;
//   //
//   //   for(unsigned i=0; i<number_of_grid_coordinates; i++)
//   //   {
//   //     Lp = npt_grid_pt->Lp_pt(i);
//   //     f = npt_grid_pt->get_box_force(i);
//   //
//   //     // update the step with the force
//   //     *Lp += step_size * time_step * f;
//   //   }
//   // }
//
//   // void O(NptGrid* npt_grid_pt)
//   // {
//   //     // dereference the box coordinates
//   //     double* Lp;
//   //     double mass = npt_grid_pt->get_mass();
//   //
//   //     // loop over the box coordinates
//   //     for(unsigned i=0; i<number_of_grid_coordinates; i++)
//   //     {
//   //         // dereference
//   //         Lp = npt_grid_pt->Lp_pt(i);
//   //
//   //         // update
//   //         *Lp = o_box_step_c * (*Lp) + o_box_step_zeta * sqrt(mass)
//   //               * gsl_ran_gaussian(box_generator[i], 1.0);
//   //     }
//   //  }
// };
// // class NPT : public Thermostat
// // {
// // private:
// //   double a_x; // box distance in a direction
// //
// //   double b_x; // box distance in b-x direction
// //   double b_y; // box distance in b-y direction
// //
// //   double xi_a_x; // momentum of box in a_x direction
// //
// //   double xi_b_x; // momentum of boc in b_x direction
// //   double xi_b_y; // mementum of box in b_y direction
// //
// //   gsl_rng* r_xi_a_x; // random seed a_x momentum
// //
// //   gsl_rng* r_xi_b_x; // random seed b_x momentum
// //   gsl_rng* r_xi_b_y; // random seed b_y momentum
// //
// //   double VOLUME; // current volume of the box
// //   double PRESSURE; // current pressur of the box
// //
// //   double MU;   // box mass
// //   double ETA;  // Langevin Dynamics Friction
// //
// //   double Time_Step; // Time_step for integration
// //
// //   FILE* testf;
// //
// //   // function which returns the pressure
// //   double get_pressure();
// //
// //   // function which returns the scaled coordinates
// //   std::vector<double> get_scaled_coord
// //   (Particle* ParPt);
// //
// //   // function which returns the scaled momentum
// //   // of the particle PartPt
// //   std::vector<double> get_scaled_momentum
// //   (Particle* PartPt);
// //
// //   // function which rescales the coordinates
// //   // of the particle back into the box. Here
// //   // A_x, B_x and B_y are the lengths of the
// //   // previous box
// //   void rescale_into_volume
// //   (Particle* PartPt,
// //    const double& A_x,
// //    const double& B_x,
// //    const double& B_y);
// //
// // public:
// //   // constructor
// //   NPT
// //   (const double& L_x,
// //    const double& L_y,
// //    const double& pressure,
// //    const double& mu,
// //    const double& eta,
// //    const double& timestep);
// //
// //   // destructor
// //   ~NPT()
// //   {
// //     // free all the random number generator memory
// //     gsl_rng_free(r_xi_a_x);
// //     gsl_rng_free(r_xi_b_x);
// //     gsl_rng_free(r_xi_b_y);
// //   };
// //
// //   // A step in volume
// //   void A
// //   (Molecule* MolPt,
// //    const double& step_size);
// //
// //   // B step in momentum
// //   void B
// //   (Molecule* MolPt,
// //    const double& step_size);
// //
// //   // O step in momentum
// //   void O
// //   (const double& kT,
// //    const double& step_size);
// //
// //   // get the volume of the box
// //   double get_volume();
// //
// //   // update the grid integrator helper
// //   void update_volume
// //   (Molecule* MolPt,
// //    Grid* GridPt);
// // };
//
// #endif
