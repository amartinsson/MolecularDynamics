#ifndef MAGNETISATIONOBSERVABLE_HPP
#define MAGNETISATIONOBSERVABLE_HPP

#include <fstream>
#include <sstream>
#include <math.h>

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "HistObservable.hpp"
#include "Molecules.hpp"

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class MagnetisationObservable : public SystemObservable
{
public:
    // empty constructor
    MagnetisationObservable(Molecule* molecule_pt, const double& mmin,
        const double& mmax, const int& N, const int& recf, const int& rect);
    ~MagnetisationObservable(){}
    // virtual function for updating
    void update();
    void update(const double& weight);
    // print function at time_index
    void print(const char* file_name,
               const unsigned& index);
    // get average
    double get_average();
    // get instant
    double get_instant();
    // observation for radial distribution
    HistObservable magDist;
private:
    // hold pointer to molecule
    Molecule* molecule_pt;
};
#endif
