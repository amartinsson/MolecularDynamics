#include "InfiniteSwitchSimulatedTempering.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
InfiniteSwitchSimulatedTempering::InfiniteSwitchSimulatedTempering(
    Molecule* molecule_pt, const double& T_min, const double& T_max,
        const unsigned& nint, const double& time_step, const double& tau)
    : InfiniteSwitch(molecule_pt, 1.0/T_max-1.0/T_max, 1.0/T_min-1.0/T_max, nint, time_step/tau)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize the force
    InfiniteSwitch::initialize_force();
}

void InfiniteSwitchSimulatedTempering::initialize_collecivet_variable()
{
    // set the collective variable as potential
    InfiniteSwitch::collective_pt = &InfiniteSwitch::molecule_pt->potential();

    // add to the map for the gradients of the potential
    for(const auto& particle : InfiniteSwitch::molecule_pt->Particles) {
        InfiniteSwitch::grad_collective_pt.insert(make_pair(particle.second,
                                                        &particle.second->f));
    }
}
