#ifndef DOUBLEINFINITESWITCHMAGTEMP_HPP
#define DOUBLEINFINITESWITCHMAGTEMP_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "DoubleInfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Tempering
 *****************************************************************************/
class DoubleInfiniteSwitchMagTemp : public DoubleInfiniteSwitch
{
public:
    DoubleInfiniteSwitchMagTemp(Molecule* molecule_pt,
        const double& T_min, const double& T_max,
            const double& b_min, const double& b_max,
                const unsigned& nTint, const unsigned& nBint,
                    const double& time_step, const double& tau);

    ~DoubleInfiniteSwitchMagTemp() {};

    // get the collective variables
    double get_collective_one();
    double get_collective_two();
    // get the gradient of the collective variables
    Vector get_collective_grad_one(Particle* particle);
    Vector get_collective_grad_two(Particle* particle);


private:
    void initialize_collecivet_variables();
    Vector dummy_gradient;

    // collective variable stuff
    const double* collective_one_pt; // pointer to collective variable
    unordered_map<Particle*, Vector*> grad_collective_one_pt;

    const double* collective_two_pt;
};

#endif
