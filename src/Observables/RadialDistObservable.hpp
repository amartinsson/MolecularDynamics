#ifndef RADIALDISTOBSERVABLE_HPP
#define RADIALDISTOBSERVABLE_HPP

#include <fstream>
#include <sstream>

#include "SystemObservable.hpp"
#include "HistObservable.hpp"

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class RadialDistObservable : public SystemObservable
{
public:
    // empty constructor
    RadialDistObservable(const double& rmin, const double& rmax, const int& N);
    // empty destructor
    ~RadialDistObservable(){};
    // virtual function for updating
    void update();
    void update(const double& r);
    // print function at time_index
    void print(const char* file_name, const double& time,
               const unsigned& index);
    // print function at time
    void print(const char* file_name, const double& time) {};

private:
    // observation for radial distribution
    HistObservable radialDist;
};
#endif
