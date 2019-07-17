#ifndef SYSTEMPOTENTIALENERGY_HPP
#define SYSTEMPOTENTIALENERGY_HPP

#include "SystemEnergy.hpp"

using namespace::std;

/******************************************************************************
                        System Potential Energy Class
 *****************************************************************************/
class SystemPotentialEnergy : public SystemEnergy
{
public:
    // constructor
    SystemPotentialEnergy(const Molecule* molecule_pt, const int& recf,
                          const int& rect);
    // destructor
    ~SystemPotentialEnergy() {/* empty */ };
    // update function
    void update();
};

#endif
