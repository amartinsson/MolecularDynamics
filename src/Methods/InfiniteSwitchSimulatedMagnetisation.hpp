#ifndef INFINITESWITCHSIMULATEDMAGNETISATION_HPP
#define INFINITESWITCHSIMULATEDMAGNETISATION_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Magnetisation
 *****************************************************************************/
class InfiniteSwitchSimulatedMagnetisation : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedMagnetisation(Molecule* molecule_pt,
        const double& M_min, const double& M_max, const unsigned& nint,
            const double& time_step, const double& tau);

private:
    void initialize_collecivet_variable();
    Vector dummy_gradient;
};

#endif
