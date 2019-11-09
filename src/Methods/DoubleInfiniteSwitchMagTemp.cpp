#include "DoubleInfiniteSwitchMagTemp.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
DoubleInfiniteSwitchMagTemp::DoubleInfiniteSwitchMagTemp(Molecule* molecule_pt,
    const double& T_min, const double& T_max,
        const double& b_min, const double& b_max,
            const unsigned& nTint, const unsigned& nBint,
                const double& time_step, const double& tau, const double& T_ref)
    : DoubleInfiniteSwitch(molecule_pt, 1.0/T_ref-1.0/T_min, 1.0/T_ref-1.0/T_max,
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
    return molecule_pt->nparticle() * (*collective_two_pt);
}

// get the gradient of the collective variables
Vector  DoubleInfiniteSwitchMagTemp::get_collective_grad_one(Particle* particle)
{
    // return the negative of the force
    return (*grad_collective_one_pt.at(particle)).neg();
}

Vector  DoubleInfiniteSwitchMagTemp::get_collective_grad_two(Particle* particle)
{
    // this is not correct, this is the derivative of the magnetisation w.r.t
    // the angle. This means that the gradient should be \nabla_{\theta_i} m
    // which is cosine: cosine of N * theta?

    return dummy_gradient * (-1.0) * sin(particle->q(0))
            * molecule_pt->nparticle();
}

double DoubleInfiniteSwitchMagTemp::get_mid_beta() {

    double beta_mid = get_mid_lambda_one();

    printf("getting mid beta = %f\n", beta_mid);

    return beta_mid;
}
