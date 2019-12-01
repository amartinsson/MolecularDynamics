#ifndef INFINITESWITCH_HPP
#define INFINITESWITCH_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"
#include "NptGrid.hpp"

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
class InfiniteSwitch
{
public:
    InfiniteSwitch(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
        const unsigned& threshold=0);
    // destructor
    ~InfiniteSwitch();
    // get the thermal rescaling factor
    void apply_force_rescaling();
    // get all the observable weights for recording in histograms
    vector<double> get_observable_weights();
    // return the number of interpolation points
    unsigned get_interpolation_points() const {return nint;}

    // must have an initialization for the collective variable
    virtual double get_collective() = 0;
    virtual Vector get_collective_grad(Particle* particle) = 0;
    virtual Matrix get_collective_grad_rot(Particle* particle) = 0;
    virtual Matrix get_collective_virial_grad() = 0;

    // initialize method to use npt or not
    void initialize_with_npt(NptGrid* npt_grid_pt);

    void print_weights(const char* filename, const double& t);

protected:
    // pointer to molecule
    Molecule* molecule_pt;

    // pointer to npt grid
    NptGrid* npt_pt;

    // initialize the force
    void initialize_force();

    // initialize with NPT
    void initialize_with_npt();

    double get_mid_lambda();

private:
    // pointer to the potential
    const double* potential_pt;

    // boolean to use npt or not
    bool with_npt;

    // vector holders
    double* gauss_weight;
    vector<double> omega_weight;
    double* lambda;

    unsigned Threshold;
    unsigned StepCount;

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
    // apply the force rescaling on all parts
    void apply_force_rescaling_all();
};

#endif
