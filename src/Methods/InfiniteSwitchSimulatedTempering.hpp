#ifndef INFINITESWITCHSIMULATEDTEMPERING_HPP
#define INFINITESWITCHSIMULATEDTEMPERING_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Tempering
 *****************************************************************************/
class InfiniteSwitchSimulatedTempering : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedTempering(Molecule* molecule_pt,
        const double& T_min, const double& T_max,
            const unsigned& nint, const double& time_step, const double& tau);

    // get the gradient of the particle
    double get_collective();
    Vector get_collective_grad(Particle* particle);
    Matrix get_collective_grad_rot(Particle* particle);
    Matrix get_collective_virial_grad();

private:
    void initialize_collecivet_variable();

    // collective variable stuff
    const double* collective_pt; // pointer to collective variable
    // this is a map between the forces and the gradient of the
    // collective variable for each particle
    unordered_map<Particle*, Vector*> grad_collective_pt;
};

#endif
