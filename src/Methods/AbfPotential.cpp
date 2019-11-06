#include "AbfPotential.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
AbfPotential::AbfPotential(
    Molecule* molecule_pt, const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step, const unsigned& threshold)
    : AdaptiveBiasingForce(molecule_pt, lambda_min, lambda_max, nint,
                                time_step, threshold)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();
}

void AbfPotential::initialize_collecivet_variable() {
    // set the collective variable as potential
    collective_pt = &AdaptiveBiasingForce::molecule_pt->potential();
}

Vector AbfPotential::get_collective_grad(Particle* particle) {
    // return the negative of the force i.e gradient
    return particle->f.neg();
}

double AbfPotential::get_collective_laplace(Particle* particle) {
    // return the laplacian of the particle
    return particle->laplace;
}

double AbfPotential::get_collective() {

    // if(BaseRef == 0.0) {
    //     BaseRef = (*collective_pt);
    //     printf("using baseref in ISST: %e\n", BaseRef);
    // }

    // return *collective_pt;
    // return *collective_pt - BaseRef;
    return *collective_pt;
}
