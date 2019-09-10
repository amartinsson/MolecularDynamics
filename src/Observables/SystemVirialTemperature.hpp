#ifndef SYSTEMVIRIALTEMPERATURE_HPP
#define SYSTEMVIRIALTEMPERATURE_HPP

#include "SystemTemperature.hpp"

using namespace::std;

/******************************************************************************
                        System Temperature Class
 *****************************************************************************/
class SystemVirialTemperature : public SystemTemperature
{
public:
    // constructor
    SystemVirialTemperature(Molecule* molecule_pt, const int& recf,
                                     const int& rect);
    // destructor
    ~SystemVirialTemperature();
    // update
    void update();
    // calculate configurational temperature
    void update_temperature();
};

#endif
