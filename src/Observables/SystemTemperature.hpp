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
    SystemTemperature(Molecule* molecule_pt);
    // destructor
    ~SystemTemperature(){};
    // set the momentum temperature
    void set_momentum_temp(){momentum = true;}
    // set the momentum temperature
    void set_config_temp(){configurational = true;}
    // instant
    double get_instant();
    // average
    double get_average();
    // update function
    void update();
    // print function at timeindex
    void print(const char* file_name, const double& time,
               const unsigned& time_index);
    // print function at time
    void print(const char* file_name, const double& time);

private:
    // local pointer to molecule object
    Molecule* system;
    // private momentum temperature
    AverageObservable* mTemp;
    // private configurational temperature
    AverageObservable* cTemp;

    // boolean for momentum temperature
    bool momentum;
    // boolean for configurational temperature
    bool configurational;

    // calculate momentum temperature
    void update_mtemp();
    // calculate configurational temperature
    void update_ctemp();
};

#endif
