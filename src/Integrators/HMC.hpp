#ifndef HMC_HPP
#define HMC_HPP

#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <unordered_map>

#include "Hamilton.hpp"
#include "MetropolisHastings.hpp"
#include "Particle.hpp"
#include "System.hpp"
#include "Array.hpp"

using namespace::std;

/******************************************************************************
                           Hamiltonain Base Class
 *****************************************************************************/
class HMC : public Hamilton, public MetropolisHastings {
public:
    HMC(const double& time_step, const double& beta,
        System* system_pt, const int& seed);
    // destructor
    ~HMC();
    // integrator
    void integrate(Molecule* molecule_pt);

protected:
    // class pointers
    // System* System_pt;
    double time_step;
    double H;
    double Hnp1;
    // system for storing positions
    unordered_map<Particle*, Vector> position;

    void log_current_position(Molecule* molecule_pt);

    void hamiltoniain_forward(Molecule* molecule_pt);

    void propsal_momentum(Particle& particle);

    double calculate_hamiltonian(const Particle& particle);
};

#endif
