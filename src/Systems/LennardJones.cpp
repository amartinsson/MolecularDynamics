#include "LennardJones.hpp"

using namespace::std;

LennardJones::LennardJones(const double& epsilon, const double& sigma)
{
    // set the force parameters
    Epsilon = epsilon;
    Sigma = sigma;
}

// destructor
LennardJones::~LennardJones()
{
    delete[] &Epsilon;
    delete[] &Sigma;
}

// compute the force
void LennardJones::compute_force(Molecule* molecule_pt)
{
    std::cout << "Error: Lennard Jones force is a pair force\n";
    exit(-1);
}

// compute the pair force
vector<double> LennardJones::compute_pair_force(Molecule* molecule_pt,
                                                Particle* particle_i,
                                                Particle* particle_j,
                                                const double& r,
                                                const double& r_x,
                                                const double& r_y)
{
    // tracker of the pair force
     vector<double> forces(2, 0.0);
     double* V = molecule_pt->potential_pt();

     // ------------ LENNARD JONES FORCE
     //
     //    F_x = - 24 * eps * (2 * sigma^12 * r^-13 - sigma^6 r^-7) * x/r;
     //    F_y = - 24 * eps * (2 * sigma^12 * r^-13 - sigma^6 r^-7) * y/r;
     //

     double dphidr = -24.* Epsilon * pow(Sigma, 6.0) *
                    (2 * pow(Sigma, 6.0) * pow(r, -14.0) - pow(r, -8.0)) * r;

   // #pragma omp atomic
   //   *PartI->f_pt(0) += 2.0 * (dphidr * r_x / r);
   //
   // #pragma omp atomic
   //   *PartI->f_pt(1) += 2.0 * (dphidr * r_y / r);
   //
   // #pragma omp atomic
   //   *PartJ->f_pt(0) -= 2.0 * (dphidr * r_x / r);
   //
   // #pragma omp atomic
   //   *PartJ->f_pt(1) -= 2.0 * (dphidr * r_y / r);

   #pragma omp atomic
     *particle_i->f_pt(0) += dphidr * r_x / r;

   #pragma omp atomic
     *particle_i->f_pt(1) += dphidr * r_y / r;

   #pragma omp atomic
     *particle_j->f_pt(0) -= dphidr * r_x / r;

   #pragma omp atomic
     *particle_j->f_pt(1) -= dphidr * r_y / r;

     // add and return the pair force
     forces[0] = dphidr * r_x / r;
     forces[1] = dphidr * r_y / r;

     // add to the potential
   #pragma omp atomic
     *V += LennardJones::get_potential(r);

     return forces;
}

// function for the Lennard Jones potential
double LennardJones::get_potential(const double& r)
{
    return 4 * Epsilon * (pow(Sigma / r, 12) - pow(Sigma / r, 6));
}
