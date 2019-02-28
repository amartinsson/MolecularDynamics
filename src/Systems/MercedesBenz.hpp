#ifndef MERCEDESBENZ_HPP
#define MERCEDESBENZ_HPP

#include "LennardJones.hpp"
#include "Matrix.hpp"

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
    vector<double> compute_pair_force(Molecule* molecule_pt,
                                      Particle* particle_i,
                                      Particle* particle_j,
                                      const double& r, const double& r_x,
                                      const double& r_y);
private:
    double Epsilon_LJ;
    double Sigma_LJ;
    double Epsilon_HB;
    double Sigma_HB;
    double R_HB;

    // Gaussian function
    double G(const double& x);
};

#endif
