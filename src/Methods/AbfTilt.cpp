#include "AbfTilt.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
AbfTilt::AbfTilt(
    Molecule* molecule_pt, const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step, const unsigned& threshold)
    : AdaptiveBiasingForce(molecule_pt, lambda_min, lambda_max, nint,
                                time_step, threshold),
        dummy_gradient(molecule_pt->dim(), 1.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();
}

void AbfTilt::initialize_collecivet_variable() {
    // HACK
    // set the collective variable as the magnetisation,
    // this is somewhat hacky and needs to be calculated in the
    // system
    // currently tested for double well for field of (x+1)
    collective_pt = &AdaptiveBiasingForce::molecule_pt->magnetisation();
}

Vector AbfTilt::get_collective_grad(Particle* particle) {
    // return the negative of the force i.e gradient
    return dummy_gradient;
}

double AbfTilt::get_collective_laplace(Particle* particle) {
    // return the laplacian of the particle
    return 0.0;
}

double AbfTilt::get_collective() {

    // if(BaseRef == 0.0) {
    //     BaseRef = (*collective_pt);
    //     printf("using baseref in ISST: %e\n", BaseRef);
    // }

    // return *collective_pt;
    // return *collective_pt - BaseRef;
    return *collective_pt;
}
