#ifndef SphereOrderObservable_HPP
#define SphereOrderObservable_HPP

#include <fstream>
#include <sstream>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "Molecules.hpp"
#include "SystemObservable.hpp"

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class SphereOrderObservable : public SystemObservable
{
public:
    // constructor
    SphereOrderObservable(Molecule* molecule_pt, const int& l_max,
                          const double& coff, const int& recf,
                            const int& rect);
    // destructor
    ~SphereOrderObservable();
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
    // holder for molecule
    Molecule* system;
    // number of spherical harmonics to calculate
    size_t lmax;
    int looplmax;
    // length of output array in sphere calculations
    int N;
    // bool to keep track of if to record step or not
    bool locRecStep;
    // gsl norm function
    gsl_sf_legendre_t norm;
    // cutoff distance for observable
    double cutoff;
    // map from particles to average observables
    unordered_map<Particle*, vector<AverageObservable*> > obsReal;
    unordered_map<Particle*, vector<AverageObservable*> > obsImag;

    // allocating all the data structure
    void allocate_unordered_maps();
    // get the polar and azimuthal polar coordinate angles
    vector<double> get_angles(const Vector& dr);
    // add observation to particle with angles
    void add_obs(Particle* particle, const vector<double>& angles);
    // calculate the order parameter for angle
    Matrix calculate_order_param(const vector<double>& angles);
    // long term observable
    AverageObservable* state;



    // map from particles to average observables
    // unordered_map<Particle*, AverageObservable*> observablemapReal;
    // unordered_map<Particle*, AverageObservable*> observablemapImag;
    // // local molecule object
    // Molecule* system;
    // // cutoff distance for observable
    // double cutoff;
    // double pr;
    // // vector with wave numbers
    // vector<Vector> k;
    // // long term observable
    // AverageObservable* state;
    //
    // // calculate the observable
    // vector<double> get_orderparam(const Vector& dr);
    // // set the wave vector
    // void set_wave();
};
#endif
