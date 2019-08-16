#ifndef CURIEWEISS_HPP
#define CURIEWEISS_HPP

#include <math.h>

#include "System.hpp"

using namespace::std;

class CurieWeiss : public System
{
public:
    // constructor
    CurieWeiss(const double& b, const unsigned& K);
    // destructor
    ~CurieWeiss();
    // compute the force
    void compute_force(Molecule* molecule_pt);
    // compute the force
    Vector compute_force(Molecule* molecule_pt, Particle* particle_i,
                Particle* particle_j, const double& r, const Vector& dr) {};
    // compute the potential
    double compute_potential(const Vector& x) {
        printf("ERROR: calling compute_potential from CurieWeiss\n");
        exit(-1);
        return 0.0;
    }
    // compute the laplacian
    double compute_laplace(const Vector& x) {
        printf("ERROR: calling compute_laplace from CurieWeiss\n");
        exit(-1);
        return 0.0;
    }

private:
    // field varialbe
    const double b;
    // number of particles variable
    const unsigned K;

    double compute_magnetisation(Molecule* molecule_pt);
    double compute_potential(const double& m);
    Vector get_force(const double& m, const Vector& theta);

    void map_periodic(Particle& particle);
};

#endif
