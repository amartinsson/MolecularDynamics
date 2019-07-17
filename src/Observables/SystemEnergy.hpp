#ifndef SYSTEMENERGY_HPP
#define SYSTEMENERGY_HPP

#include "AverageObservable.hpp"
#include "Molecules.hpp"
#include "SystemObservable.hpp"

using namespace::std;

/******************************************************************************
                        System Energy Class
 *****************************************************************************/
class SystemEnergy : public SystemObservable
{
public:
    // constructor
    SystemEnergy(const Molecule* molecule_pt, const int& recf,
                 const int& rect);
    // destructor
    ~SystemEnergy();
    // instant
    double get_instant();
    // average
    double get_average();
    // update function
    void update();

protected:
    // private momentum temperature
    AverageObservable* Energy;
    // Pointer to molecule object
    const Molecule* molecule_pt;

    // calculate kinetic energy
    double kinetic();
};

#endif
