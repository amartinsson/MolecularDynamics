#ifndef ORDEROBSERVABLE_HPP
#define ORDEROBSERVABLE_HPP

#include <fstream>
#include <sstream>
#include <math.h>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "Molecules.hpp"
#include "SystemObservable.hpp"

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class OrderObservable : public SystemObservable
{
public:
    // constructor
    OrderObservable(Molecule* molecule_pt, const double& part_rad);
    // destructor
    ~OrderObservable();
    // clear the observations
    void clear();
    // update function
    void update();
    void update(Particle* particle_1, Particle* particle_2,
                const double& r, const Vector& dr);
    // print function at time_index
    void print(const char* file_name, const double& time,
               const unsigned& index);
    // print function at time
    void print(const char* file_name, const double& time);
    // get the instant order value
    double get_instant();
    // get the average order value
    double get_average();


private:
    // map from particles to average observables
    unordered_map<Particle*, AverageObservable*> observablemapReal;
    unordered_map<Particle*, AverageObservable*> observablemapImag;
    // local molecule object
    Molecule* system;
    // cutoff distance for observable
    double cutoff;
    double pr;
    // vector with wave numbers
    vector<Vector> k;
    // long term observable
    AverageObservable* state;

    // calculate the observable
    vector<double> get_orderparam(const Vector& dr);
    // set the wave vector
    void set_wave();
};
#endif
