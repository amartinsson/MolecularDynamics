#ifndef RADIALDISTOBSERVABLE_HPP
#define RADIALDISTOBSERVABLE_HPP


#include <fstream>
#include <sstream>
#include <math.h>

#include "AverageObservable.hpp"
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
    RadialDistObservable(const double& rmin, const double& rmax, const int& N,
                         const int& recf, const int& rect,
                         const int& nparticles);
    // empty destructor
    ~RadialDistObservable(){delete number_density;}
    // bump the counter for RecStep
    void bump_recstep(const double& V);
    // virtual function for updating
    void update();
    void update(const double& r);
    // print function at time_index
    void print(const char* file_name, const double& time,
               const unsigned& index);
    // get average
    double get_average();
    // get instant
    double get_instant();
    // reset
    void reset(const int& recf, const int& rect) {

        SystemObservable::reset(recf, rect);

        number_density->clear();

        radialDist.clear();
    }

private:
    // observation for radial distribution
    HistObservable radialDist;
    // local recording step
    bool localRecStep;
    // local number of particles
    int number_of_particles;
    // aveerage number density
    AverageObservable* number_density;

    // get the radial distribution
    double get_rpdf_3d(const int& i);
};
#endif
