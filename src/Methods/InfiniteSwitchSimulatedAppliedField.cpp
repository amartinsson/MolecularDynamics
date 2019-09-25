#include "InfiniteSwitchSimulatedAppliedField.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated AppliedField
//  *****************************************************************************/
InfiniteSwitchSimulatedAppliedField::InfiniteSwitchSimulatedAppliedField(
    Molecule* molecule_pt, const double& b_min, const double& b_max,
        const unsigned& nint, const double& time_step, const double& tau)
    : InfiniteSwitch(molecule_pt, b_min, b_max, nint, time_step/tau),
        dummy_gradient(molecule_pt->dim(), 1.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize the force
    InfiniteSwitch::initialize_force();
}

void InfiniteSwitchSimulatedAppliedField::initialize_collecivet_variable()
{
    // set the collective variable as the magnetisation,
    // this is somewhat hacky and needs to be calculated in the
    // system
    // currently tested for double well for field of (x+1)
    collective_pt = &InfiniteSwitch::molecule_pt->magnetisation();
}

Vector InfiniteSwitchSimulatedAppliedField::get_collective_grad(Particle* particle)
{
    return dummy_gradient;
}

double InfiniteSwitchSimulatedAppliedField::get_collective() {
    return *collective_pt;
}
