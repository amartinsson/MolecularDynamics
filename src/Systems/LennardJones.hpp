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

    Vector compute_force(Molecule* molecule_pt,
                         Particle* particle_i,
                         Particle* particle_j,
                         const double& r,
                         const Vector& dr);
     // calculate the potential
     double compute_potential(const Vector& x) {/* empty */};
     double compute_laplace(const Vector& x) {/* empty */};
     
protected:
    // force parameters
    double Epsilon;
    double Sigma;

    // function for the Lennard Jones potential
    double get_potential(const double& r);
    // function which returns the force scalar
    double get_force_scalar(const double& r);
};

#endif
