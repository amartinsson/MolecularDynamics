#ifndef SphereOrderObservable_HPP
#define SphereOrderObservable_HPP

#include <fstream>
#include <sstream>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "Molecules.hpp"
#include "SystemObservable.hpp"

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
@article{auer_quantitative_2004,
	title = {{QUANTITATIVE} {PREDICTION} {OF} {CRYSTAL}-{NUCLEATION} {RATES} {FOR} {SPHERICAL} {COLLOIDS}: {A} {Computational} {Approach}},
	volume = {55},
	shorttitle = {{QUANTITATIVE} {PREDICTION} {OF} {CRYSTAL}-{NUCLEATION} {RATES} {FOR} {SPHERICAL} {COLLOIDS}},
	url = {https://doi.org/10.1146/annurev.physchem.55.091602.094402},
	doi = {10.1146/annurev.physchem.55.091602.094402},
	abstract = {This review discusses the recent progress that has been made in the application of computer simulations to study crystal nucleation in colloidal systems. We discuss the concept and the numerical methods that allow for a quantitative prediction of crystal-nucleation rates. The computed nucleation rates are predicted from first principles and can be directly compared with experiments. These techniques have been applied to study crystal nucleation in hard-sphere colloids, polydisperse hard-sphere colloids, weakly charged or slightly soft colloids, and hard-sphere colloids that are confined between two-plane hard walls.},
	number = {1},
	urldate = {2019-11-10},
	journal = {Annual Review of Physical Chemistry},
	author = {Auer, Stefan and Frenkel, Daan},
	year = {2004},
	pmid = {15117256},
	pages = {333--361},
	file = {Full Text PDF:/home/anton/Zotero/storage/DFRDZDTA/Auer and Frenkel - 2004 - QUANTITATIVE PREDICTION OF CRYSTAL-NUCLEATION RATE.pdf:application/pdf}
}
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

    void print_traj(const char* file_name, const unsigned& index);

private:
    // holder for molecule
    Molecule* system;
    // number of spherical harmonics to calculate
    unsigned lmax;
    int looplmax;
    // bool to keep track of if to record step or not
    bool locRecStep;
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
};
#endif
