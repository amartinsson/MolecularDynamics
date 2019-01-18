#ifndef LENNARDJONES_HPP
#define LENNARDJONES_HPP

#include <math.h>

#include "System.hpp"

using namespace::std;

class LennardJones : public System
{
public:
    LennardJones(const double& epsilon, const double& sigma);
    // destructor
    ~LennardJones();
    // compute the force
    void compute_force(Molecule* molecule_pt);
    // compute the pair force
    vector<double> compute_pair_force(Molecule* molecule_pt,
                                      Particle* particle_i,
                                      Particle* particle_j,
                                      const double& r, const double& r_x,
                                      const double& r_y);
protected:
    // force parameters
    double Epsilon;
    double Sigma;

    // function for the Lennard Jones potential
    double get_potential(const double& r);
};

#endif
