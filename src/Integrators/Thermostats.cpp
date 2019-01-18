// // --------------------------------------------------------------------------------------------------------------------//
// //                                              THERMOSTAT CLASS
// // --------------------------------------------------------------------------------------------------------------------//
//
// #include "Thermostats.hpp"
//
// // --------------------------------------------------------------------------------------------------------------------//
// // Particle A step
// void NptThermostat::A_two(Particle* particle_pt, NptGrid* npt_grid_pt,
//                           const double& step_size)
// {
//     //---------------------------------------------------- PARTICLE INTEGRATION
//     // holder for particle positon
//     double* q = NULL;
//
//     // dereference the box size
//     double lx1 = *npt_grid_pt->L_pt(0);
//     double ly1 = *npt_grid_pt->L_pt(1);
//     double ly2 = *npt_grid_pt->L_pt(2);
//
//     // get the momentum of the particle
//     double px = *particle_pt->p_pt(0);
//     double py = *particle_pt->p_pt(1);
//
//     // get the momentum of the particle
//     double mx = *particle_pt->m_pt(0);
//     double my = *particle_pt->m_pt(1);
//
//     // loop over the number of dimensions
//     for(unsigned j=0; j<dimension; j++)
//     {
//         // dereference the particle position
//         q = particle_pt->q_pt(j);
//
//         if(j==0)
//         {
//             (*q) += step_size * time_step *((1.0 - ly1*ly1/(lx1*ly2)) * px/mx
//                                             + ly1/ly2 * py/my);
//         }
//         else if(j==1)
//         {
//             (*q) += step_size * time_step *(py/my - ly1/lx1 * px/mx);
//         }
//     }
//
//     // -------------------------------------------------------- BOX INTEGRATION
//     if(box_integration == true)
//     {
//         // update the accumulated momentum in the box
//         npt_grid_pt->update_accumulted_momentum();
//         double acc_mom = 0.0;
//         double* Lp = NULL;
//
//         // loop over all the box coordinates
//         for(unsigned i=0; i<number_of_grid_coordinates;i++)
//         {
//             acc_mom = npt_grid_pt->get_accumulated_momentum(i);
//             Lp = npt_grid_pt->Lp_pt(i);
//
//             if(i != 2)
//                 (*Lp) += step_size * time_step * (acc_mom / lx1);
//             else
//             {
//                 double acc_mom_xy = npt_grid_pt->get_accumulated_momentum(1);
//                 (*Lp) += step_size * time_step *
//                          (acc_mom / ly2 - ly1/(ly2*lx1) * acc_mom_xy);
//             }
//         }
//
//         // double m0 = npt_grid_pt->get_accumulated_momentum(0);
//         // double m1 = npt_grid_pt->get_accumulated_momentum(1);
//         // double m2 = npt_grid_pt->get_accumulated_momentum(2);
//         //
//         // printf("Accummulated momentum %2.3f %2.3f %2.3f\n", m0, m1, m2);
//
//         // set boolean to false
//         box_integration = false;
//     }
// }
//
// // Particle B step
// void NptThermostat::B(Particle* particle_pt, NptGrid* npt_grid_pt,
//                       const double& step_size)
// {
//     //---------------------------------------------------- PARTICLE INTEGRATION
//     // holder for particle momentum
//     double* p = NULL;
//
//     // dereference the box size
//     double lx1 = *npt_grid_pt->L_pt(0);
//     double ly1 = *npt_grid_pt->L_pt(1);
//     double ly2 = *npt_grid_pt->L_pt(2);
//
//     // get the forces on the particle
//     double fx = *particle_pt->f_pt(0);
//     double fy = *particle_pt->f_pt(1);
//
//     for(unsigned j=0; j<dimension; j++)
//     {
//         // derefernce the particle momentum
//         p = particle_pt->p_pt(j);
//
//         if(j==0)
//         {
//             (*p) += step_size * time_step * ((1.0 - ly1*ly1/(lx1*ly2)) * fx - ly1/lx1 * fy);
//         }
//         else if(j==1)
//         {
//             (*p) += step_size * time_step * (fy + ly1/ly2 * fx);
//         }
//     }
//
//     // -------------------------------------------------------- BOX INTEGRATION
//     if(box_integration == true)
//     {
//         double* Lp = NULL;
//         double virial = 0.0;
//
//         for(unsigned i=0; i<number_of_grid_coordinates; i++)
//         {
//             Lp = npt_grid_pt->Lp_pt(i);
//             virial = npt_grid_pt->get_virial(i);
//
//             if(i == 0)
//                 *Lp -= step_size * time_step * (virial + target_pressure * ly2);
//             else if(i == 1)
//                 *Lp -= step_size * time_step * virial;
//             else if(i == 2)
//                 *Lp -= step_size * time_step * (virial + target_pressure * lx1);
//         }
//
//         // set boolean to false
//         box_integration = false;
//     }
// }
//
// // Particle O step
// void NptThermostat::O(Particle* particle_pt, NptGrid* npt_grid_pt)
// {
//     //---------------------------------------------------- PARTICLE INTEGRATION
//     // particlce mass
//     double  m  = 0.0;
//
//     // holder for particle momentum
//     double* p = NULL;
//
//     // loop over all the dimensions
//     for(unsigned j=0; j<dimension; j++)
//     {
//         m = *particle_pt->m_pt(j);
//         p = particle_pt->p_pt(j);
//
//         (*p) = o_part_step_c * (*p)
//                      + o_part_step_zeta * sqrt(m)
//                        * gsl_ran_gaussian(part_generator[j], 1.0);
//     }
//
//     // -------------------------------------------------------- BOX INTEGRATION
//     if(box_integration == true)
//     {
//         // dereference the box coordinates
//         double* Lp;
//         double box_mass = npt_grid_pt->get_mass();
//
//         // loop over the box coordinates
//         for(unsigned i=0; i<number_of_grid_coordinates; i++)
//         {
//             // dereference
//             Lp = npt_grid_pt->Lp_pt(i);
//
//             // update
//             *Lp = o_box_step_c * (*Lp) + o_box_step_zeta * sqrt(box_mass)
//                   * gsl_ran_gaussian(box_generator[i], 1.0);
//         }
//
//         // set boolean to false
//         box_integration = false;
//     }
// }
