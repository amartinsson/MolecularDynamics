#ifndef ADAPTIVEBIASINGFORCE_HPP
#define ADAPTIVEBIASINGFORCE_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "LegendreRuleFast.hpp"
#include "InfiniteSwitch.hpp"
#include "Molecules.hpp"
#include "NptGrid.hpp"

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
class AdaptiveBiasingForce
{
public:
    AdaptiveBiasingForce(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
        const unsigned& threshold);
    // destructor
    ~AdaptiveBiasingForce();
    // get the thermal rescaling factor
    void apply_force_rescaling();
    // get all the observable weights for recording in histograms
    vector<double> get_observable_weights();
    // return the number of interpolation points
    unsigned get_interpolation_points() const {return nint;}

    // must have an initialization for the collective variable
    virtual double get_collective() = 0;
    virtual Vector get_collective_grad(Particle* particle) = 0;
    virtual double get_collective_laplace(Particle* particle) = 0;

    // initialize method to use npt or not
    void initialize_with_npt(NptGrid* npt_grid_pt);

protected:
    // pointer to molecule
    Molecule* molecule_pt;

    // pointer to npt grid
    NptGrid* npt_pt;

    // initialize the force
    void initialize_force();

private:
    // pointer to the potential
    const double* potential_pt;

    unsigned Threshold;
    unsigned StepCount;

    // vector average holders
    AverageObservable* partition_estimate;
    AverageObservable* force_estimate;
    HistObservable ilocator;

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
    // void update_hull();
    // update the adaptive force
    void update_adaptive_force(const int& index);
    // update the weights
    // void update_weights();
    // must have private calculate force rescaling
    // double calculate_lambda_bar();
    // rescaler the weights directly
    // void normalize_weights();
    // apply the force rescaling on all parts
    void apply_force_rescaling(const int& index);
};

#endif
