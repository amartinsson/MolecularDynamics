#ifndef BAB_HPP
#define BAB_HPP

#include "Hamilton.hpp"

using namespace::std;

/******************************************************************************
                           BAB Hamiltonain Integrator
 *****************************************************************************/
class BAB : public Hamilton {
public:
    BAB(const double& time_step, System* system_pt);
    // destructor
    ~BAB() {};
    // integrator
    void integrate(Molecule* molecule_pt);
private:
    double time_step;
};

#endif
