#include "SimulatedAnnealing.hpp"

using namespace::std;
/******************************************************************************
                               Simulated Tempering
 *****************************************************************************/
SimulatedAnnealing::SimulatedAnnealing(const double& tmax,
    const double& crate, const unsigned& nsteps, Molecule* molecule_pt)
        : tmax(tmax), crate(crate), nsteps(nsteps),
            molecule_pt(molecule_pt) {}

// update the temperature according to MH
void SimulatedAnnealing::update_temperature(const unsigned& step)
{
    double fraction = (double)step / (double)nsteps;

    if(fraction < 1.0) {
        if(crate == 0) {
            lin_cooling(fraction);
        }
        else {
            exp_cooling(fraction);
        }
    }
    else {
        printf("WARNING: Calling SimTemp::update_temperature with 1\n");
    }

    molecule_pt->set_beta(beta);
}

// return the current temperature of the simulation
double SimulatedAnnealing::get_beta()
{
    return beta;
}

// exponential cooling rate
void SimulatedAnnealing::exp_cooling(const double& f)
{
    double T = tmax * (exp(-crate * f) - exp(-crate)) / (1.0 - exp(-crate));

    beta = 1.0 / T;
}

// linear cooling rate
void SimulatedAnnealing::lin_cooling(const double& f)
{
    double T = tmax * (1.0 - f);

    beta = 1.0 / T;
}
