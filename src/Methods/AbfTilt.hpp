#ifndef ABFTILT_HPP
#define ABFTILT_HPP

#include "AdaptiveBiasingForce.hpp"

using namespace::std;

/******************************************************************************
                Adaptive biasing Force with Tilt Simulation
 *****************************************************************************/
class AbfTilt : public AdaptiveBiasingForce
{
public:
    AbfTilt(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
        const unsigned& threshold);
    // destructor
    ~AbfTilt();

    // must have an initialization for the collective variable
    double get_collective();
    Vector get_collective_grad(Particle* particle);
    double get_collective_laplace(Particle* particle);

private:
    // must have an initialization for the collective variable
    void initialize_collecivet_variable();
    // collective variable stuff
    const double* collective_pt;
    // collective variable hessian
    const double* collective_laplace_pt;

    // collective dummy_gradient
    Vector dummy_gradient;
};

#endif
