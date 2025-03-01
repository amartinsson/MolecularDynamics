#ifndef INFINITESWITCHSIMULATEDMAGNETISATION_HPP
#define INFINITESWITCHSIMULATEDMAGNETISATION_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Magnetisation
 *****************************************************************************/
class InfiniteSwitchSimulatedMagnetisation : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedMagnetisation(Molecule* molecule_pt,
        const double& b_min, const double& b_max, const unsigned& nint,
            const double& time_step, const double& tau);

    // get the gradient of the particle
    double get_collective();

    // get the gradient of the collective variable for particle
    Vector get_collective_grad(Particle* particle);

private:
    // initializes the collective variable to point to the correct position
    void initialize_collecivet_variable();

    // collective variable stuff
    const double* collective_pt; // pointer to collective variable
    
    // store a dummy gradient of ones
    Vector dummy_gradient;
};

#endif
