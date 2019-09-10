#ifndef INFINITESWITCHSIMULATEDTEMPERING_HPP
#define INFINITESWITCHSIMULATEDTEMPERING_HPP

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
class InfiniteSwitchSimulatedTempering : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedTempering(Molecule* molecule_pt,
        const double& T_min, const double& T_max,
            const unsigned& nint, const double& time_step, const double& tau);

private:
    void initialize_collecivet_variable();
};

#endif
