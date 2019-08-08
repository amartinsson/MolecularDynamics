#ifndef SIMULATEDTEMPERING_HPP
#define SIMULATEDTEMPERING_HPP

#include <vector>
#include <omp.h>

#include "Generator.hpp"
#include "AverageObservable.hpp"
#include "Molecules.hpp"

using namespace::std;
/******************************************************************************
                               Simulated Tempering
 *****************************************************************************/
class SimulatedTempering
{
public:
    SimulatedTempering(const double& tmin, const double& tmax,
                       const double& n_temperatures,
                       const unsigned& mod_switch,
                       const int& seed);
    // destructor
    ~SimulatedTempering();
    // print method info to stream
    void get_method_info();
    // update the temperature according to MH
    void update_temperature(Molecule* molecule_pt, const unsigned& step);
    // return the current temperature index
    unsigned get_temperature_index();
    // return the current temperature of the simulation
    double get_temperature();
    // get the i^th beta
    double get_beta(const unsigned& i);
    // set the temperature weight for i^th component
    void set_beta_weight(const double& weight, const unsigned& i);

protected:
    // inverse temperature limits
    double beta_min;
    double beta_max;
    // current temperature index
    unsigned current_index;
    unsigned number_of_temperatures;
    // hold the switching frequency of the algorithm
    unsigned switch_frequency;
    // holder for all the temperatures
    std::vector<double> beta;
    std::vector<double> Weights;
    std::vector<bool> state_visited;
    // keeper of average potential energy
    std::vector<AverageObservable*> ave_V;
    // uniform random variable
    //gsl_rng* uniform;
    UniformGenerator uniform_gen;

    // rescaling the weights such that they sum to one
    void normalise_weights();
    // rescaling of the momentum
    void rescale_momentum(Molecule* molecule_pt, const int& shift);
};

#endif
