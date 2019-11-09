#ifndef DOUBLEINFINITESWITCH_HPP
#define DOUBLEINFINITESWITCH_HPP

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
class DoubleInfiniteSwitch
{
public:
    DoubleInfiniteSwitch(Molecule* molecule_pt,
        const double& lambda_one_min, const double& lambda_one_max,
        const double& lambda_two_min, const double& lambda_two_max,
        const unsigned& nOneint, const unsigned& nTwoint,
        const double& time_step);
    // destructor
    ~DoubleInfiniteSwitch();
    // get the thermal rescaling factor
    void apply_force_rescaling();
    // get all the observable weights for recording in histograms
    Matrix get_observable_weights();
    // get the collective variables
    virtual double get_collective_one() = 0;
    virtual double get_collective_two() = 0;
    // get the gradient of the collective variables
    virtual Vector get_collective_grad_one(Particle* particle) = 0;
    virtual Vector get_collective_grad_two(Particle* particle) = 0;

    // return the number of interpolation points
    vector<unsigned> get_interpolation_points() const {

        vector<unsigned> retv = {nOneint, nTwoint};

        return retv;
    }


protected:
    // pointer to molecule
    // const Molecule* molecule_pt;
    Molecule* molecule_pt;

    // initialize the force
    void initialize_force();

    double get_mid_lambda_one();
    double get_mid_lambda_two();

private:
    // pointer to the potential
    const double* potential_pt;

    // vector holders
    double* gauss_one_weight;
    double* lambda_one;

    // vector holders
    double* gauss_two_weight;
    double* lambda_two;

    // matrix of omega weights
    Matrix omega_weight;

    // matrix lookup for hull_estimate
    unordered_map<int, vector<AverageObservable*> > hull_estimate;

    // double holders
    const double lambda_one_min;
    const double lambda_one_max;

    const double lambda_two_min;
    const double lambda_two_max;

    const double time_step;
    const double beta;

    // unsigned holders
    const unsigned nOneint;
    const unsigned nTwoint;

    // must have an initialization for the collective variable
    virtual void initialize_collecivet_variables() = 0;

    // update the hull estimate
    void update_hull();
    // update the weights
    void update_weights();
    // must have private calculate force rescaling
    Vector calculate_lambda_bar();
    // rescaler the weights directly
    void normalize_weights();
};

#endif
