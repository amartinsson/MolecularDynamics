#include "InfiniteSwitchSimulatedMagnetisation.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Magnetisation
//  *****************************************************************************/
InfiniteSwitchSimulatedMagnetisation::InfiniteSwitchSimulatedMagnetisation(
    Molecule* molecule_pt, const double& M_min, const double& M_max,
        const unsigned& nint, const double& time_step, const double& tau)
    : InfiniteSwitch(molecule_pt, 1.0/M_max-1.0/M_max, 1.0/M_min-1.0/M_max,
        nint, time_step/tau), dummy_gradient(molecule_pt->dim(), 1.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize the force
    InfiniteSwitch::initialize_force();
}

void InfiniteSwitchSimulatedMagnetisation::initialize_collecivet_variable()
{
    // set the collective variable as potential
    InfiniteSwitch::collective_pt =
        &InfiniteSwitch::molecule_pt->magnetisation();

    // add to the map for the gradients of the potential
    for(const auto& particle : InfiniteSwitch::molecule_pt->Particles) {
        InfiniteSwitch::grad_collective_pt.insert(make_pair(particle.second,
                                                        &dummy_gradient));
    }
}
