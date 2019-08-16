#ifndef INFINITESWITCHSIMULATEDTEMPERING_HPP
#define INFINITESWITCHSIMULATEDTEMPERING_HPP

#include <vector>
#include <omp.h>

#include "AverageObservable.hpp"
#include "InfiniteSwitch.hpp"
#include "LegendreRuleFast.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                    Infinite Switch Simulated Tempering
 *****************************************************************************/
class InfiniteSwitchSimulatedTempering : public InfiniteSwitch
{
public:
    InfiniteSwitchSimulatedTempering(Molecule* molecule_pt,
        const double& lambda_min, const double& lambda_max,
            const unsigned& nint, const double& time_step, const double& tau);

private:
    void initialize_collecivet_variable();
};

// /******************************************************************************
//                     Infinite Switch Simulated Tempering
//  *****************************************************************************/
// class InfiniteSwitchSimulatedTempering
// {
// public:
//     InfiniteSwitchSimulatedTempering(Molecule* molecule_pt, const double& tmin,
//                                      const double& tmax, const unsigned& nint,
//                                      const double& t_step, const double& tau);
//     // destructor
//     ~InfiniteSwitchSimulatedTempering();
//     // get the thermal rescaling factor
//     void apply_force_rescaling(Molecule* molecule_pt);
//     // get the i^th interpolation point
//     double get_beta(const unsigned& i);
//     // get all the observable weights for recording in histograms
//     vector<double> get_observable_weights(Molecule* molecule_pt);
//     // set all the beta weights
//     void set_beta_weight(const vector<double>& weight, Molecule* molecule_pt);
//
// private:
//     // vector holders
//     double* beta;
//     double* gauss_weight;
//     double* beta_weight;
//
//     // vector average holders
//     AverageObservable* partition_estimate;
//
//     // double holders
//     double beta_min;
//     double beta_max;
//     double time_step;
//     double thermal_force_scaling;
//     double time_step_scaling;
//     double potential_scale;
//
//     // unsigned holders
//     unsigned number_of_interpolation_points;
//
//     // calculate the force rescaling
//     void calculate_force_rescaling(Molecule* molecule_pt);
//     // initialise the force
//     void init_force(Molecule* molecule_pt);
//     // learn the beta weights
//     void learn_beta_weights(Molecule* molecule_pt);
//     // rescaler the weights directly
//     void normalise_weights();
//     // update the thermal force resacling
//     void update_force_rescaling(Molecule* molecule_pt);
// };

#endif
