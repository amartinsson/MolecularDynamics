#include "PressureInfiniteSwitchSimulation.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
PressureInfiniteSwitchSimulation::PressureInfiniteSwitchSimulation(
    Molecule* molecule_pt, NptGrid* npt_pt, const double& P_min,
        const double& P_max, const unsigned& nint, const double& time_step,
            const double& tau)
    : InfiniteSwitch(molecule_pt,
        molecule_pt->beta() * (npt_pt->get_target_pressure() - P_max),
            molecule_pt->beta() * (npt_pt->get_target_pressure() - P_min),
                nint, time_step/tau), dummy_grad(molecule_pt->dim(), 0.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize with npt
    InfiniteSwitch::initialize_with_npt(npt_pt);

    // initialize the force
    InfiniteSwitch::initialize_force();
}

void PressureInfiniteSwitchSimulation::initialize_collecivet_variable()
{
    // just return the volume from the npt_pt
}

Vector PressureInfiniteSwitchSimulation
    ::get_collective_grad(Particle* particle)
{
    return dummy_grad;
}

double PressureInfiniteSwitchSimulation::get_collective()
{
    return npt_pt->S.det();
}

Matrix PressureInfiniteSwitchSimulation::get_collective_virial_grad()
{
    return npt_pt->virial.neg();
}
