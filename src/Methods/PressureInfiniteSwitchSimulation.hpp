#ifndef PRESSUREINFINITESWITCHSIMULATION_HPP
#define PRESSUREINFINITESWITCHSIMULATION_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Tempering
 *****************************************************************************/
class PressureInfiniteSwitchSimulation : public InfiniteSwitch
{
public:
    PressureInfiniteSwitchSimulation(Molecule* molecule_pt, NptGrid* npt_pt,
        const double& P_min, const double& P_max, const unsigned& nint,
            const double& time_step, const double& tau);

    // get the gradient of the particle
    double get_collective();
    Vector get_collective_grad(Particle* particle);
    Matrix get_collective_grad_rot(Particle* particle);
    Matrix get_collective_virial_grad();

private:
    void initialize_collecivet_variable();

    // collective variable stuff
    const double* collective_pt;
    // this should be a dim length zero vector
    Vector dummy_grad;
    Matrix dummy_grad_rot;
};

#endif
