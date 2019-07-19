#ifndef HAMILTON_HPP
#define HAMILTON_HPP
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

#include "Integrator.hpp"
#include "Particle.hpp"
#include "System.hpp"
#include "Array.hpp"

using namespace::std;

/******************************************************************************
                           Hamiltonain Base Class
 *****************************************************************************/
class Hamilton : public Integrator {
public:
    Hamilton(System* system_pt);
    // destructor
    ~Hamilton();
    // integrator
    void integrate(Molecule* molecule_pt);

protected:
    // class pointers
    System* System_pt;
    
    void A(Particle& particle, const double& h);
    void B(Particle& particle, const double& h);
};

#endif
