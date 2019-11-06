#ifndef ABFPOTENTIAL_HPP
#define ABFPOTENTIAL_HPP

#include "AdaptiveBiasingForce.hpp"

using namespace::std;

/******************************************************************************
                Adaptive biasing Force with Potential Simulation
 *****************************************************************************/
class AbfPotential : public AdaptiveBiasingForce
{
public:
    AbfPotential(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
        const unsigned& threshold);
    // destructor
    ~AbfPotential();

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
};

#endif
