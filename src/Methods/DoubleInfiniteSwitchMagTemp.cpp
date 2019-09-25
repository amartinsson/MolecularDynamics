#include "DoubleInfiniteSwitchMagTemp.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
DoubleInfiniteSwitchMagTemp::DoubleInfiniteSwitchMagTemp(Molecule* molecule_pt,
    const double& T_min, const double& T_max,
        const double& b_min, const double& b_max,
            const unsigned& nTint, const unsigned& nBint,
                const double& time_step, const double& tau)
    : DoubleInfiniteSwitch(molecule_pt, 1.0/T_min-1.0/T_min, 1.0/T_min-1.0/T_max,
            b_min, b_max, nTint, nBint, time_step/tau),
                dummy_gradient(molecule_pt->dim(), 1.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variables();

    // initialize the force
    DoubleInfiniteSwitch::initialize_force();
}

void DoubleInfiniteSwitchMagTemp::initialize_collecivet_variables()
{
    // set the collective variable as potential
    collective_one_pt = &molecule_pt->potential();

    collective_two_pt = &molecule_pt->magnetisation();

    // add to the map for the gradients of the potential
    for(const auto& particle : DoubleInfiniteSwitch::molecule_pt->Particles) {
        grad_collective_one_pt.insert(make_pair(
                particle.second, &particle.second->f));
    }
}

// get the collective variables
double  DoubleInfiniteSwitchMagTemp::get_collective_one()
{
    return *collective_one_pt;
}

double  DoubleInfiniteSwitchMagTemp::get_collective_two()
{
    return *collective_two_pt;
}

// get the gradient of the collective variables
Vector  DoubleInfiniteSwitchMagTemp::get_collective_grad_one(Particle* particle)
{
    // return the negative of the force
    return (*grad_collective_one_pt.at(particle)).neg();
}

Vector  DoubleInfiniteSwitchMagTemp::get_collective_grad_two(Particle* particle)
{
    return dummy_gradient;
}
