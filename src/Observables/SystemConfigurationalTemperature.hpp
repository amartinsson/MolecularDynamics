#ifndef SYSTEMCONFIGURATIONALTEMPERATURE_HPP
#define SYSTEMCONFIGURATIONALTEMPERATURE_HPP

#include "SystemTemperature.hpp"

using namespace::std;

/******************************************************************************
                        System Temperature Class
 *****************************************************************************/
class SystemConfigurationalTemperature : public SystemTemperature
{
public:
    // constructor
    SystemConfigurationalTemperature(Molecule* molecule_pt, const int& recf,
                                     const int& rect);
    // destructor
    ~SystemConfigurationalTemperature();
    // get instant
    double get_instant();
    // get average 
    double get_average();
    // calculate configurational temperature
    void update_temp();

private:
    AverageObservable* nablaSquare;
    AverageObservable* laplace;
};

#endif
