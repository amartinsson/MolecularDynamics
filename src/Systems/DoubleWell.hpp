#ifndef DOUBLEWELL_HPP
#define DOUBLEWELL_HPP

#include <math.h>

#include "System.hpp"

using namespace::std;

class DoubleWell : public System
{
public:
    DoubleWell(const double& a, const double& b);
    // destructor
    ~DoubleWell() {/* empty */};
    // compute the force
    void compute_force(Molecule* molecule_pt);
    // compute the force
    Vector compute_force(Molecule* molecule_pt,
                         Particle* particle_i,
                         Particle* particle_j,
                         const double& r,
                         const Vector& dr);
    double compute_potential(const Vector& x) {return get_potential(x);}
    double compute_laplace(const Vector& x) {return get_laplace(x);}
protected:
    // force parameters
    double a;
    double b;

    // function for the potential
    double get_potential(const double& x);
    double get_potential(const Vector& x);
    // function for the laplace
    double get_laplace(const double& x);
    double get_laplace(const Vector& x);
    // function which returns the force scalar
    double get_force(const double& x);
    Vector get_force(const Vector& x);
};

#endif
