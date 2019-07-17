#ifndef MERCEDESBENZ_HPP
#define MERCEDESBENZ_HPP

#include "LennardJones.hpp"
#include "Array.hpp"
#include <math.h>

using namespace::std;

/******************************************************************************
                              Mercedes Benz
    Mercedes Benz system is an implemntation of a system described in:
    https://arxiv.org/pdf/1304.3232.pdf
 *****************************************************************************/
class MercedesBenz : public LennardJones
{
public:
    MercedesBenz(const double& epsilon_LJ, const double& sigma_LJ,
                 const double& epsilon_HB, const double& sigma_HB,
                 const double& r_HB);
    //destructor
    ~MercedesBenz();
    // compute the force
    void compute_force(Molecule* molecule_pt);
    // compute the pair force
    Vector compute_force(Molecule* molecule_pt,
                         Particle* particle_i,
                         Particle* particle_j,
                         const double& r,
                         const Vector& dr);

private:
    double Epsilon_HB;
    double Sigma_HB;
    double R_HB;

    // Gaussian function
    double G(const double& x);
};

#endif
