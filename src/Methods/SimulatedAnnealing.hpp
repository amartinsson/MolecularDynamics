#ifndef SIMULATEDANNEALING_HPP
#define SIMULATEDANNEALING_HPP

#include <vector>
#include <omp.h>
#include <math.h>

#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                               Simulated Tempering
 *****************************************************************************/
class SimulatedAnnealing
{
public:
    SimulatedAnnealing(const double& tmax, const double& crate,
        const unsigned& nsteps, Molecule* molecule_pt);
    // destructor
    ~SimulatedAnnealing();
    // print method info to stream
    void get_method_info();
    // update the temperature according to MH
    void update_temperature(const unsigned& step);
    // return the current temperature of the simulation
    double get_beta();

protected:
    // inverse temperature limits
    double tmax;
    unsigned nsteps;
    // current temperture
    double beta;
    double crate;

    Molecule* molecule_pt;

    void exp_cooling(const double& fraction);
    void lin_cooling(const double& fraction);
};

#endif
