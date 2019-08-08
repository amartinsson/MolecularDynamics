#ifndef ONGULATEDTEMPERING_HPP
#define ONGULATEDTEMPERING_HPP

#include "SimulatedTempering.hpp"

using namespace::std;
/******************************************************************************
                               Ongulated Tempering
 *****************************************************************************/
class OngulatedTempering : public SimulatedTempering {
public:
    // constructor
    OngulatedTempering(const double& tmin, const double& tmax,
                       const double& n_temperatures,
                       const unsigned& mod_switch,
                       const int& seed);
    // destructor
    ~OngulatedTempering(){};
    // update the temperature
    void update_temperature(Molecule* molecule_pt, const unsigned& step);

};

#endif
