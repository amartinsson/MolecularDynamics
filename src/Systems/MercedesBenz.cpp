#include "MercedesBenz.hpp"

using namespace::std;

/******************************************************************************
                              Mercedes Benz
    Mercedes Benz system is an implemntation of a system described in:
    https://arxiv.org/pdf/1304.3232.pdf
 *****************************************************************************/
MercedesBenz::MercedesBenz(const double& epsilon_LJ, const double& sigma_LJ,
                           const double& epsilon_HB, const double& sigma_HB,
                           const double& r_HB)
                           : LennardJones(epsilon_LJ, sigma_LJ)
{
    // set all the force parameters
    Epsilon_LJ = epsilon_LJ;
    Sigma_LJ = sigma_LJ;
    Epsilon_HB = epsilon_HB;
    Sigma_HB = sigma_HB;
    R_HB = R_HB;
}

MercedesBenz::~MercedesBenz()
{
    delete &Epsilon_LJ;
    delete &Sigma_LJ;
    delete &Epsilon_HB;
    delete &Sigma_HB;
    delete &R_HB;
}

// compute the force
void MercedesBenz::compute_force(Molecule* molecule_pt)
{
    std::cout << "Error: Mercedes Benz force is a pair force\n";
    exit(-1);
}

// compute the pair force
vector<double> MercedesBenz::compute_pair_force(Molecule* molecule_pt,
                                                Particle* particle_i,
                                                Particle* particle_j,
                                                const double& r,
                                                const double& r_x,
                                                const double& r_y)
{
      // holders for computing the forces
      // between the particle_icles
      std::vector<double> lj_forces(2, 0.0);
      std::vector<double> forces(2, 0.0);

      // Lennard Jones Force and Potential
      lj_forces = LennardJones::compute_pair_force(molecule_pt, particle_i,
                                                   particle_j, r, r_x, r_y);

      forces[0] += lj_forces[0];
      forces[1] += lj_forces[1];

    #pragma omp atomic
      *molecule_pt->potential_pt() += LennardJones::get_potential(r);

      // ------------------------------ HYDROGEN BOND ----------------------------------- //

      // Helpers and derivatives with respect to x and y
      std::vector<double> h_i(3, 0.0);
      std::vector<double> h_j(3, 0.0);

      std::vector<double> dh_i_dx_i(3,0.0);
      std::vector<double> dh_j_dx_i(3,0.0);

      std::vector<double> dG_dh_i(3,0.0);
      std::vector<double> dG_dh_j(3,0.0);

      std::vector<double> dh_i_dy_i(3,0.0);
      std::vector<double> dh_j_dy_i(3,0.0);

      std::vector<double> dh_i_dQ_00(3,0.0);
      std::vector<double> dh_i_dQ_01(3,0.0);

      std::vector<double> dh_i_dQ_10(3,0.0);
      std::vector<double> dh_i_dQ_11(3,0.0);

      std::vector<double> dh_j_dQ_00(3,0.0);
      std::vector<double> dh_j_dQ_01(3,0.0);

      std::vector<double> dh_j_dQ_10(3,0.0);
      std::vector<double> dh_j_dQ_11(3,0.0);

      double dphiHB_dx = 0.0;
      double dphiHB_dy = 0.0;

      double dphiHB_dQ_i_00 = 0.0;
      double dphiHB_dQ_i_01 = 0.0;

      double dphiHB_dQ_i_10 = 0.0;
      double dphiHB_dQ_i_11 = 0.0;

      double dphiHB_dQ_j_00 = 0.0;
      double dphiHB_dQ_j_01 = 0.0;

      double dphiHB_dQ_j_10 = 0.0;
      double dphiHB_dQ_j_11 = 0.0;

      // derivatives with respect to r
      double dG_dr = -(r-R_HB) / (Sigma_HB*Sigma_HB) * G(r-R_HB);
      double dr_dx = -r_x/r;
      double dr_dy = -r_y/r;

      // unit vector in r direction
      std::vector<double> u(2,0.0);
      u[0] = r_x/r;
      u[1] = r_y/r;

      // local storage for particle_icle
      double Q_i[2][2];
      double Q_j[2][2];

      double arm_i[2][3];
      double arm_j[2][3];

      // read in orientation from memory
      // for(unsigned k=0; k<2; k++)
      //   for(unsigned l=0;l<2;l++)
      //   {
      //       // particle_icle i
      //       Q_i[k][l] = *particle_i->Q_pt(k, l);
      //
      //       // particle_icle j
      //       Q_j[k][l] = *particle_j->Q_pt(k, l);
      //   }
      for(unsigned k=0; k<3; k++)
      {
            for(unsigned l=0;l<2;l++)
            {
                if(k<2)
                {
                    // particle_icle i
                    //Q_i[k][l] = *particle_i->Q_pt(k, l);
                    Q_i[k][l] = particle_i->Q(k, l);

                    // particle_icle j
                //    Q_j[k][l] = *particle_j->Q_pt(k, l);
                    Q_j[k][l] = particle_j->Q(k, l);
                }

                // read in particle_icle i
                arm_i[l][k] = particle_i->arm(l, k);

                // read in particle_icle j
                arm_j[l][k] = particle_j->arm(l, k);
            }
        }

      // // read in arms from memory
      // for(unsigned k=0; k<3; k++)
      // {
      //     // read in particle_icle i
      //     arm_i[0][k] = particle_i->arm(0, k);
      //     arm_i[1][k] = particle_i->arm(1, k);
      //
      //     // read in particle_icle j
      //     arm_j[0][k] = particle_j->arm(0, k);
      //     arm_j[1][k] = particle_j->arm(1, k);
      // }

      for(unsigned k=0; k<3; k++)
      {
          h_i[k] = (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) * u[0]
                 + (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) * u[1];

          h_j[k] = (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) * u[0]
               + (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) * u[1];

          dh_i_dx_i[k] = r_x * h_i[k] / (r * r)
                       - (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) / r;

          dh_j_dx_i[k] = r_x * h_j[k] / (r * r)
                       - (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) / r;

          dh_i_dy_i[k] = r_y * h_i[k] / (r * r)
                       - (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) / r;

          dh_j_dy_i[k] = r_y * h_j[k] / (r * r)
                       - (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) / r;

          dh_i_dQ_00[k] = u[0] *  arm_i[0][k];
          dh_i_dQ_01[k] = u[1] *  arm_i[0][k];

          dh_i_dQ_10[k] = u[0] *  arm_i[1][k];
          dh_i_dQ_11[k] = u[1] *  arm_i[1][k];

           // matrix elements for the Torque particle_icle J
          dh_j_dQ_00[k] = u[0] *  arm_j[0][k];
          dh_j_dQ_01[k] = u[1] *  arm_j[0][k];

          dh_j_dQ_10[k] = u[0] *  arm_j[1][k];
          dh_j_dQ_11[k] = u[1] *  arm_j[1][k];

      }

      // Helpers for exponential calculation
      std::vector<double> G_hi_m1 = {G(h_i[0]-1.0), G(h_i[1]-1.0), G(h_i[2]-1.0)};
      std::vector<double> G_hj_p1 = {G(h_j[0]+1.0), G(h_j[1]+1.0), G(h_j[2]+1.0)};

      // optimiser helpers
      double sum_G_hj_p1 = 0.0;
      double G_r_rHB = G(r-R_HB);

      double sum_dphiHB_dx = 0.0;
      double sum_dphiHB_dy = 0.0;

      double sum_dphiHB_dQ_j_00 = 0.0;
      double sum_dphiHB_dQ_j_01 = 0.0;
      double sum_dphiHB_dQ_j_10 = 0.0;
      double sum_dphiHB_dQ_j_11 = 0.0;

      // Assign the values from above
      for(unsigned k=0;k<3;k++)
        {
            // calculate sum
            sum_G_hj_p1 += G_hj_p1[k];

          // h_i
        //   h_i[k] =
    	// ((*particle_i->Q_pt(0,0)) * particle_i->arm(0,k) + (*particle_i->Q_pt(1,0)) * particle_i->arm(1,k)) * u[0] +
    	// ((*particle_i->Q_pt(0,1)) * particle_i->arm(0,k) + (*particle_i->Q_pt(1,1)) * particle_i->arm(1,k)) * u[1];
          // h_i[k] = (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) * u[0]
          //        + (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) * u[1];

           // h_j
         //   h_j[k] =
    	 // ((*particle_j->Q_pt(0,0)) * particle_j->arm(0,k) + (*particle_j->Q_pt(1,0)) * particle_j->arm(1,k)) * u[0] +
    	 // ((*particle_j->Q_pt(0,1)) * particle_j->arm(0,k) + (*particle_j->Q_pt(1,1)) * particle_j->arm(1,k)) * u[1];
           // h_j[k] = (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) * u[0]
    	   //        + (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) * u[1];

         //   dG_dh_i[k] =
    	 // -(h_i[k] - 1.0) / (Sigma_HB * Sigma_HB) * G(h_i[k] - 1.0);
         //
         //   // dG_dh_j
         //   dG_dh_j[k] =
    	 // -(h_j[k] + 1.0) / (Sigma_HB * Sigma_HB) * G(h_j[k] + 1.0);
           dG_dh_i[k] = -(h_i[k] - 1.0) / (Sigma_HB * Sigma_HB) * G_hi_m1[k];

           // dG_dh_j
           dG_dh_j[k] = -(h_j[k] + 1.0) / (Sigma_HB * Sigma_HB) * G_hj_p1[k];


         //   // dh_i_dx_i
         // //   dh_i_dx_i[k] =
    	 // // r_x * h_i[k] / (r * r)
    	 // // - ((*particle_i->Q_pt(0,0)) * particle_i->arm(0,k) + (*particle_i->Q_pt(1,0)) * particle_i->arm(1,k)) / r;
         //   dh_i_dx_i[k] = r_x * h_i[k] / (r * r)
         //                - (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) / r;
         //
         //   // dh_j_dx_i
         // //   dh_j_dx_i[k] =
    	 // // r_x * h_j[k] / (r * r)
    	 // // - ((*particle_j->Q_pt(0,0)) * particle_j->arm(0,k) + (*particle_j->Q_pt(1,0)) * particle_j->arm(1,k)) / r;
         //    dh_j_dx_i[k] = r_x * h_j[k] / (r * r)
         //                 - (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) / r;
         //
         //   // dh_i_dy_i
         // //   dh_i_dy_i[k] =
    	 // // r_y * h_i[k] / (r * r)
    	 // // - ((*particle_i->Q_pt(0,1)) * particle_i->arm(0,k) + (*particle_i->Q_pt(1,1)) * particle_i->arm(1,k)) / r;
         //    dh_i_dy_i[k] = r_y * h_i[k] / (r * r)
         //                 - (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) / r;
         //
         //   // dh_j_dy_i
         // //   dh_j_dy_i[k] =
    	 // // r_y * h_j[k] / (r * r)
    	 // // - ((*particle_j->Q_pt(0,1)) * particle_j->arm(0,k) + (*particle_j->Q_pt(1,1)) * particle_j->arm(1,k)) / r;
         //    dh_j_dy_i[k] = r_y * h_j[k] / (r * r)
         //                 - (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) / r;

           // matrix elements for the Torque particle_icle I
           // dh_i_dQ_00[k] = u[0] *  particle_i->arm(0,k);
           // dh_i_dQ_01[k] = u[1] *  particle_i->arm(0,k);
           //
           // dh_i_dQ_10[k] = u[0] *  particle_i->arm(1,k);
           // dh_i_dQ_11[k] = u[1] *  particle_i->arm(1,k);
           //
           // // matrix elements for the Torque particle_icle J
           // dh_j_dQ_00[k] = u[0] *  particle_j->arm(0,k);
           // dh_j_dQ_01[k] = u[1] *  particle_j->arm(0,k);
           //
           // dh_j_dQ_10[k] = u[0] *  particle_j->arm(1,k);
           // dh_j_dQ_11[k] = u[1] *  particle_j->arm(1,k);
           // dh_i_dQ_00[k] = u[0] *  arm_i[0][k];
           // dh_i_dQ_01[k] = u[1] *  arm_i[0][k];
           //
           // dh_i_dQ_10[k] = u[0] *  arm_i[1][k];
           // dh_i_dQ_11[k] = u[1] *  arm_i[1][k];
           //
           // // matrix elements for the Torque particle_icle J
           // dh_j_dQ_00[k] = u[0] *  arm_j[0][k];
           // dh_j_dQ_01[k] = u[1] *  arm_j[0][k];
           //
           // dh_j_dQ_10[k] = u[0] *  arm_j[1][k];
           // dh_j_dQ_11[k] = u[1] *  arm_j[1][k];

           // calculate sum
           sum_dphiHB_dx += dG_dh_j[k] * dh_j_dx_i[k];
           sum_dphiHB_dy += dG_dh_j[k] * dh_j_dy_i[k];

           sum_dphiHB_dQ_j_00 += dG_dh_j[k] * dh_j_dQ_00[k];
           sum_dphiHB_dQ_j_01 += dG_dh_j[k] * dh_j_dQ_01[k];
           sum_dphiHB_dQ_j_10 += dG_dh_j[k] * dh_j_dQ_10[k];
           sum_dphiHB_dQ_j_11 += dG_dh_j[k] * dh_j_dQ_11[k];
         }

       // k over the arms
      //#pragma omp simd
       for(unsigned k=0;k<3;k++)
         {
             // ---------------- HYDROGEN BOND LINEAR --------------- //
             dphiHB_dx += sum_G_hj_p1 * (dG_dr * dr_dx * G_hi_m1[k]
                                        + G_r_rHB * dG_dh_i[k] * dh_i_dx_i[k])
                        + sum_dphiHB_dx * G_r_rHB * G_hi_m1[k];

             dphiHB_dy += sum_G_hj_p1 * (dG_dr * dr_dy * G_hi_m1[k]
                                        + G_r_rHB * dG_dh_i[k] * dh_i_dy_i[k])
                        + sum_dphiHB_dy * G_r_rHB * G_hi_m1[k];

            // --------------- HYDROGEN BOND TORQUE I -------------- //
            dphiHB_dQ_i_00 += dG_dh_i[k] * dh_i_dQ_00[k] * sum_G_hj_p1;
            dphiHB_dQ_i_01 += dG_dh_i[k] * dh_i_dQ_01[k] * sum_G_hj_p1;
            dphiHB_dQ_i_10 += dG_dh_i[k] * dh_i_dQ_10[k] * sum_G_hj_p1;
            dphiHB_dQ_i_11 += dG_dh_i[k] * dh_i_dQ_11[k] * sum_G_hj_p1;

            // --------------- HYDROGEN BOND TORQUE J -------------- //
            dphiHB_dQ_j_00 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_00;
            dphiHB_dQ_j_01 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_01;
            dphiHB_dQ_j_10 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_10;
            dphiHB_dQ_j_11 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_11;

            // ---------------- HYDROGEN POTENTIAL --------------- //
        #pragma omp atomic
            *molecule_pt->potential_pt() += Epsilon_HB * G_r_rHB * G_hi_m1[k] * sum_G_hj_p1;

         //   for(unsigned l=0;l<3;l++)
    	 // {
    	 //   // ---------------- HYDROGEN BOND LINEAR --------------- //
    	 //   // dphiHB_dx +=
    	 //   //     Epsilon_HB * dG_dr     * dr_dx         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
    	 //   //   + Epsilon_HB * G(r-R_HB) * dG_dh_i[k]    * dh_i_dx_i[k]  * G(h_j[l]+1.0)
    	 //   //   + Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dx_i[l];
         //   //
    	 //   // dphiHB_dy +=
    	 //   //     Epsilon_HB * dG_dr     * dr_dy         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
    	 //   //   + Epsilon_HB * G(r-R_HB) * dG_dh_i[k]    * dh_i_dy_i[k]  * G(h_j[l]+1.0)
    	 //   //   + Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dy_i[l];
         //   // dphiHB_dx +=
         //   //     dG_dr     * dr_dx         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
         //   //   + G(r-R_HB) * dG_dh_i[k]    * dh_i_dx_i[k]  * G(h_j[l]+1.0)
         //   //   + G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dx_i[l];
         //   //
         //   // dphiHB_dy +=
         //   //      dG_dr     * dr_dy         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
         //   //   + G(r-R_HB) * dG_dh_i[k]    * dh_i_dy_i[k]  * G(h_j[l]+1.0)
         //   //   +  G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dy_i[l];
         //
         //     // dphiHB_dx += dG_dr * dr_dx * G_hi_m1[k] * G_hj_p1[l]
         //     //   +  G_r_rHB * dG_dh_i[k]    * dh_i_dx_i[k]  * G_hj_p1[l]
         //     //   +  G_r_rHB * G_hi_m1[k] * dG_dh_j[l]    * dh_j_dx_i[l];
         //     //
         //     // dphiHB_dy +=
         //     //      dG_dr     * dr_dy         * G_hi_m1[k]  * G_hj_p1[l]
         //     //   +  G_r_rHB * dG_dh_i[k]    * dh_i_dy_i[k]  * G_hj_p1[l]
         //     //   +  G_r_rHB *G_hi_m1[k]  * dG_dh_j[l]    * dh_j_dy_i[l];
         //
         //   // dphiHB_dx += G_r_rHB * G_hi_m1[k] * dG_dh_j[l] * dh_j_dx_i[l];
         //   //
         //   // dphiHB_dy += G_r_rHB * G_hi_m1[k] * dG_dh_j[l] * dh_j_dy_i[l];
         //
    	 //   // --------------- HYDROGEN BOND TORQUE I -------------- //
    	 //   // dphiHB_dQ_i_00 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * dG_dh_i[k] * dh_i_dQ_00[k] * G(h_j[l]+1.0);
         //   //
    	 //   // dphiHB_dQ_i_01 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * dG_dh_i[k] * dh_i_dQ_01[k] * G(h_j[l]+1.0);
         //   //
    	 //   // dphiHB_dQ_i_10 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * dG_dh_i[k] * dh_i_dQ_10[k] * G(h_j[l]+1.0);
         //   //
    	 //   // dphiHB_dQ_i_11 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * dG_dh_i[k] * dh_i_dQ_11[k] * G(h_j[l]+1.0);
         //   // dphiHB_dQ_i_00 += dG_dh_i[k] * dh_i_dQ_00[k] * G(h_j[l]+1.0);
         //   //
         //   // dphiHB_dQ_i_01 += dG_dh_i[k] * dh_i_dQ_01[k] * G(h_j[l]+1.0);
         //   //
         //   // dphiHB_dQ_i_10 += dG_dh_i[k] * dh_i_dQ_10[k] * G(h_j[l]+1.0);
         //   //
         //   // dphiHB_dQ_i_11 += dG_dh_i[k] * dh_i_dQ_11[k] * G(h_j[l]+1.0);
         //
    	 //   // --------------- HYDROGEN BOND TORQUE J -------------- //
    	 //   // dphiHB_dQ_j_00 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_00[l];
         //   //
    	 //   // dphiHB_dQ_j_01 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_01[l];
         //   //
    	 //   // dphiHB_dQ_j_10 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_10[l];
         //   //
    	 //   // dphiHB_dQ_j_11 +=
    	 //   //   Epsilon_HB * G(r-R_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_11[l];
         //   // dphiHB_dQ_j_00 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_00[l];
         //   //
    	 //   // dphiHB_dQ_j_01 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_01[l];
         //   //
    	 //   // dphiHB_dQ_j_10 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_10[l];
         //   //
    	 //   // dphiHB_dQ_j_11 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_11[l];
         //
         //   // dphiHB_dQ_j_00 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_00[l];
         //   // dphiHB_dQ_j_01 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_01[l];
         //   // dphiHB_dQ_j_10 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_10[l];
         //   // dphiHB_dQ_j_11 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_11[l];
         //
    	 //   // // ---------------- HYDROGEN POTENTIAL --------------- //
    	 //   // *molecule_pt->potential_pt() += Epsilon_HB * G(r - R_HB) * G(h_i[k]-1) * G(h_j[l]+1);
    	 // }
         }

    //    // add force for particle_icle I
    // #pragma omp atomic
    //     *particle_i->f_pt(0) -= dphiHB_dx;
    //
    // #pragma omp atomic
    //     *particle_i->f_pt(1) -= dphiHB_dy;
    //
    //    // add force for particle_icle J
    // #pragma omp atomic
    //     *particle_j->f_pt(0) += dphiHB_dx;
    //
    // #pragma omp atomic
    //     *particle_j->f_pt(1) += dphiHB_dy;
    //
    // // add up the forces between particle_icles
    // forces[0] -= dphiHB_dx;
    // forces[1] -= dphiHB_dy;

       // add force for particle_icle I
    #pragma omp atomic
        *particle_i->f_pt(0) -= Epsilon_HB * dphiHB_dx;

    #pragma omp atomic
        *particle_i->f_pt(1) -= Epsilon_HB * dphiHB_dy;

       // add force for particle_icle J
    #pragma omp atomic
        *particle_j->f_pt(0) += Epsilon_HB * dphiHB_dx;

    #pragma omp atomic
        *particle_j->f_pt(1) += Epsilon_HB * dphiHB_dy;

    // add up the forces between particle_icles
    forces[0] -= Epsilon_HB * dphiHB_dx;
    forces[1] -= Epsilon_HB * dphiHB_dy;


       // add torque to overall torque
    // #pragma omp atomic
    //    *particle_i->tau_pt(0) +=
    //      *particle_i->Q_pt(0,1) * dphiHB_dQ_i_00 + *particle_i->Q_pt(1,1) * dphiHB_dQ_i_10
    //      - *particle_i->Q_pt(0,0) * dphiHB_dQ_i_01 - *particle_i->Q_pt(1,0) * dphiHB_dQ_i_11;
    //
    // #pragma omp atomic
    //    *particle_j->tau_pt(0) +=
    //      *particle_j->Q_pt(0,1) * dphiHB_dQ_j_00 + *particle_j->Q_pt(1,1) * dphiHB_dQ_j_10
    //      - *particle_j->Q_pt(0,0) * dphiHB_dQ_j_01 - *particle_j->Q_pt(1,0) * dphiHB_dQ_j_11;
    #pragma omp atomic
       *particle_i->tau_pt(0) += Epsilon_HB * G_r_rHB * ( Q_i[0][1] * dphiHB_dQ_i_00
                                                   + Q_i[1][1] * dphiHB_dQ_i_10
                                                   - Q_i[0][0] * dphiHB_dQ_i_01
                                                   - Q_i[1][0] * dphiHB_dQ_i_11);

    #pragma omp atomic
       *particle_j->tau_pt(0) += Epsilon_HB * G_r_rHB * ( Q_j[0][1] * dphiHB_dQ_j_00
                                                   + Q_j[1][1] * dphiHB_dQ_j_10
                                                   - Q_j[0][0] * dphiHB_dQ_j_01
                                                   - Q_j[1][0] * dphiHB_dQ_j_11);


         // deleta arrays
         // delete [] Q_i;
         // delete [] Q_j;
         // delete [] arm_i;
         // delete [] arm_j;

        // return the forces between the particle_icles
        return forces;
}

// Gaussian function
double MercedesBenz::G(const double& x)
{
    return exp(-0.5 * (x * x) / (Sigma_HB * Sigma_HB));
}
