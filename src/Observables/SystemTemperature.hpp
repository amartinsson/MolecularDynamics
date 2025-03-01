#ifndef SYSTEMTEMPERATURE_HPP
#define SYSTEMTEMPERATURE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "Molecules.hpp"

using namespace::std;

/******************************************************************************
                        System Temperature Class
 *****************************************************************************/
class SystemTemperature : public SystemObservable
{
public:
    // constructor
    SystemTemperature(Molecule* molecule_pt, const int& recf, const int& rect);
    // destructor
    ~SystemTemperature(){};
    // instant
    double get_instant();
    // average
    double get_average();
    // update function
    void update();
    // update the temperature -- can be overridden by other temperatures
    void update_temperature();

    void reset(const int& recf, const int& rect) {

        SystemObservable::reset(recf, rect);

        Temp->clear();
    }

protected:
    // local pointer to molecule object
    Molecule* system;
    // private momentum temperature
    AverageObservable* Temp;
};

#endif
