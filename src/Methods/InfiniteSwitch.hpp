#ifndef INFINITESWITCH_HPP
#define INFINITESWITCH_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
class InfiniteSwitch
{
public:
    InfiniteSwitch(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step);
    // destructor
    ~InfiniteSwitch();
    // get the thermal rescaling factor
    void apply_force_rescaling();
    // get all the observable weights for recording in histograms
    vector<double> get_observable_weights();
    // return the number of interpolation points
    unsigned get_interpolation_points() const {return nint;}

protected:
    // pointer to molecule
    // const Molecule* molecule_pt;
    Molecule* molecule_pt;
    // collective variable stuff
    const double* collective_pt; // pointer to collective variable
    unordered_map<Particle*, Vector*> grad_collective_pt;

    // initialize the force
    void initialize_force();

private:
    // pointer to the potential
    const double* potential_pt;

    // vector holders
    double* gauss_weight;
    vector<double> omega_weight;
    double* lambda;

    // vector average holders
    AverageObservable* hull_estimate;

    // double holders
    const double lambda_min;
    const double lambda_max;
    const double time_step;
    const double beta;

    // unsigned holders
    const unsigned nint;

    // must have an initialization for the collective variable
    virtual void initialize_collecivet_variable() = 0;

    // update the hull estimate
    void update_hull();
    // update the weights
    void update_weights();
    // must have private calculate force rescaling
    double calculate_lambda_bar();
    // rescaler the weights directly
    void normalize_weights();
};

#endif
