// #ifndef _METHODS_HPP_INCLUDED
// #define _METHODS_HPP_INCLUDED
//
// #include <vector>
// #include <omp.h>
// #include <gsl/gsl_integration.h>
// #include <numeric>
// #include <algorithm>
//
// #include "Molecules.hpp"
// #include "Observables.hpp"
//
// class Methods
// {
// public:
//   // empty constructor for method
//   Methods(){};
//
//   // emtpy destructor for methods
//   ~Methods(){};
//
//   // get force rescaling must be implementd by all methods
//   virtual void get_method_info() = 0;
// };
//
// // ---------------------------------------------------------------------------//
// //                     Infinite Swap Simulated Tempering
// // ---------------------------------------------------------------------------//
// class ISST : public Methods
// {
// public:
//     // constructor for the ISST method
//     ISST(Molecule* MolPt, const double& tmin, const double& tmax,
//          const double& nInt, const double& t_step, const double& alpha);
//
//     // destructor for the ISST method
//     ~ISST();
//
//     // dummy method
//     void get_method_info();
//
//     // get all the observable weights for recording in histograms
//     std::vector<double> get_observable_weights(Molecule* MolPt);
//
//     // get the i^th interpolation point
//     double get_beta(const unsigned& i);
//
//     // set all the beta weights
//     void set_beta_weight(const std::vector<double>& weight,
//                          Molecule* molecule_pt);
//
//     // update the thermal force resacling
//     void update_force_rescaling(Molecule* MolPt);
//
//     // get the thermal rescaling factor
//     void apply_force_rescaling(Molecule* molecule_pt);
//
//     // rescaler the weights directly
//     void normalise_weights();
//
// private:
//
//   // vector holders
//   std::vector<double> beta;
//   std::vector<double> gauss_weight;
//   std::vector<double> beta_weight;
//
//   // vector average holders
//   std::vector<AverageObservable*> partition_estimate;
//
//   // double holders
//   double beta_min;
//   double beta_max;
//
//   double time_step;
//   double thermal_force_scaling;
//
//   double time_step_scaling;
//
//   // unsigned holders
//   unsigned number_of_interpolation_points;
//
//   // initialise the force rescaling
//   void init_force_rescaling(Molecule* MolPt);
//
//   // initialise the force
//   void init_force(Molecule* MolPt);
//
//   // calculate the force rescaling
//   void calculate_force_rescaling(Molecule* molecule_pt);
//
//   // learn the beta weights
//   void learn_beta_weights(Molecule* molecule_pt);
//
// public:
//
//
// };
//
// // ---------------------------------------------------------------------------//
// //                              Simulated Tempering
// // ---------------------------------------------------------------------------//
// class SimulatedTempering : public Methods
// {
// public:
//     SimulatedTempering(const double& tmin, const double& tmax,
//                        const double& n_temperatures,
//                        const unsigned& mod_switch);
//
//     ~SimulatedTempering(){ gsl_rng_free(uniform); }
//
//     // print method info to stream
//     void get_method_info();
//
//     // update the temperature according to MH
//     void update_temperature(Molecule* molecule_pt, const unsigned& step);
//
//     // return the current temperature index
//     unsigned get_temperature_index();
//
//     // return the current temperature of the simulation
//     double get_temperature();
//
//     // get the i^th beta
//     double get_beta(const unsigned& i);
//
//     // set the temperature weight for i^th component
//     void set_beta_weight(const double& weight, const unsigned& i);
//
// private:
//     // inverse temperature limits
//     double beta_min;
//     double beta_max;
//
//     // current temperature index
//     unsigned current_index;
//     unsigned number_of_temperatures;
//
//     // hold the switching frequency of the algorithm
//     unsigned switch_frequency;
//
//     // holder for all the temperatures
//     std::vector<double> beta;
//     std::vector<double> Weights;
//     std::vector<bool> state_visited;
//
//     // keeper of average potential energy
//     std::vector<AverageObservable*> ave_V;
//
//     // uniform random variable
//     gsl_rng* uniform;
//
//     // rescaling the weights such that they sum to one
//     void normalise_weights();
//
//     // rescaling of the momentum
//     void rescale_momentum(Molecule* molecule_pt, const int& shift);
// };
//
//
// #endif
