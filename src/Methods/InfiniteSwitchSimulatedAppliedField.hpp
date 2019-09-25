#ifndef INFINITESWITCHSIMULATEDAPPLIEDFIELD_HPP
#define INFINITESWITCHSIMULATEDAPPLIEDFIELD_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Applied Field
 *****************************************************************************/
class InfiniteSwitchSimulatedAppliedField : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedAppliedField(Molecule* molecule_pt,
        const double& b_min, const double& b_max, const unsigned& nint,
            const double& time_step, const double& tau);

    // get the gradient of the particle
    double get_collective();
    
    // get the gradient of the collective variable for particle
    Vector get_collective_grad(Particle* particle);

private:
    void initialize_collecivet_variable();

    // collective variable stuff
    const double* collective_pt; // pointer to collective variable

    // dummy gradient of just ones
    Vector dummy_gradient;
};

#endif
