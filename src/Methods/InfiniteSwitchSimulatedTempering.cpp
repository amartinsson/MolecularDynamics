#include "InfiniteSwitchSimulatedTempering.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
InfiniteSwitchSimulatedTempering::InfiniteSwitchSimulatedTempering(
    Molecule* molecule_pt, const double& T_min, const double& T_max,
        const unsigned& nint, const double& time_step, const double& tau)
    : InfiniteSwitch(molecule_pt, 1.0/T_min-1.0/T_min, 1.0/T_min-1.0/T_max,
            nint, time_step/tau), BaseRef(0.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize the force
    // InfiniteSwitch::initialize_force();

    // BaseRef = get_collective();
    // printf("using baseref in ISST: %e\n", BaseRef);
}

void InfiniteSwitchSimulatedTempering::initialize_collecivet_variable()
{
    // set the collective variable as potential
    collective_pt = &InfiniteSwitch::molecule_pt->potential();
}

Vector InfiniteSwitchSimulatedTempering
    ::get_collective_grad(Particle* particle) {

    // return the negative of the force i.e gradient
    return particle->f.neg();
}

Matrix InfiniteSwitchSimulatedTempering
    ::get_collective_grad_rot(Particle* particle) {

    // return the negative of tau i.e gradient
    return particle->tau.neg();
}

double InfiniteSwitchSimulatedTempering::get_collective() {

    if(BaseRef == 0.0) {
        BaseRef = (*collective_pt);
        printf("using baseref in ISST: %e\n", BaseRef);
    }

    // return *collective_pt;
    return *collective_pt - BaseRef;
}

Matrix InfiniteSwitchSimulatedTempering::get_collective_virial_grad() {
    return npt_pt->virial.neg();
}
