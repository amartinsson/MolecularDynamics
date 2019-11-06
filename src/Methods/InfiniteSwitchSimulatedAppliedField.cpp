#include "InfiniteSwitchSimulatedAppliedField.hpp"

using namespace::std;

// /******************************************************************************
//                     Infinite Switch Simulated AppliedField
//  *****************************************************************************/
InfiniteSwitchSimulatedAppliedField::InfiniteSwitchSimulatedAppliedField(
    Molecule* molecule_pt, const double& b_min, const double& b_max,
        const unsigned& nint, const double& time_step, const double& tau)
    : InfiniteSwitch(molecule_pt, b_min, b_max, nint, time_step/tau),
        dummy_gradient(molecule_pt->dim(), 1.0)
{
    // initialize the collective variables
    this->initialize_collecivet_variable();

    // initialize the force
    InfiniteSwitch::initialize_force();
}

void InfiniteSwitchSimulatedAppliedField::initialize_collecivet_variable()
{
    // set the collective variable as the magnetisation,
    // this is somewhat hacky and needs to be calculated in the
    // system
    // currently tested for double well for field of (x+1)
    collective_pt = &InfiniteSwitch::molecule_pt->magnetisation();
}

Vector InfiniteSwitchSimulatedAppliedField::get_collective_grad(Particle* particle)
{
    return dummy_gradient;
}

double InfiniteSwitchSimulatedAppliedField::get_collective() {
    return *collective_pt;
}
//
// // get the thermal rescaling factor
// void InfiniteSwitchSimulatedAppliedField::apply_force_rescaling()
// {
//     // update the weights at current position
//     update_weights();
//
//     if(StepCount > Threshold) {
//         // apply the rescaling
//         apply_force_rescaling_all();
//     }
//     else {
//         StepCount++;
//
//     }
// }
//
// // get the thermal rescaling factor
// void InfiniteSwitchSimulatedAppliedField::apply_force_rescaling_all()
// {
//     // rescaling
//     double lambda_bar = calculate_lambda_bar();
//
//     // printf("lambda_bar = %f\n", lambda_bar);
//     double kT = 1.0 / beta;
//
//     for(auto& particle : molecule_pt->Particles) {
//         // get the gradient of the collective variable
//         Vector grad = get_collective_grad(particle.second);
//
//         // add to the rescaling of the force
//         particle.second->f += grad * kT * lambda_bar;
//
//         if(particle.second->rigid_body() == true) {
//             particle.second->tau += get_collective_grad_rot(particle.second) *
//                                         kT * lambda_bar;
//         }
//     }
//
//     if(with_npt) {
//         // get the dimension
//         unsigned dim = molecule_pt->dim();
//
//         // get collective virial grad
//         Matrix grad = get_collective_virial_grad();
//
//         npt_pt->virial += grad * kT * lambda_bar;
//
//         // // rescale the virial function
//         // for(unsigned i=0; i<dim; i++) {
//         //     for(unsigned j=i; j<dim; j++) {
//         //         npt_pt->virial(i, j) += grad(i, j) * kT * lambda_bar;
//         //     }
//         // }
//     }
// }
//
// // must have private calculate force rescaling
// double InfiniteSwitchSimulatedAppliedField::calculate_lambda_bar()
// {
//     // do integrals
//     double lambda_bar_lower = 0.0;
//     double lambda_bar_upper = 0.0;
//
//     double V = *potential_pt;
//     //
//     // // add to the potential if we're using npt
//     // if(with_npt) {
//     //     // V += *pressure_pt * (*cell_pt).det();
//     //     V += npt_pt->get_target_pressure() * npt_pt->S.det();
//     // }
//
//     double theta = get_collective();
//
//     for(unsigned i=0; i<nint; i++) {
//         // version one
//         // // either of these will work, first should increase stability
//         // // of exp calculation
//         double measure = omega_weight[i] * exp(-beta * V + lambda[i] * theta);
//         // double measure = omega_weight[i] * exp(lambda[i] * theta);
//
//         lambda_bar_upper += gauss_weight[i] * lambda[i] * measure;
//         lambda_bar_lower += gauss_weight[i] * measure;
//     }
//
//     return lambda_bar_upper / lambda_bar_lower;
// }
//
// // update the hull estimate
// void InfiniteSwitchSimulatedAppliedField::update_hull()
// {
//     // dereference values
//     double theta = get_collective();
//
//     // double integral = 0.0;
//     //
//     // for(unsigned j=0; j<nint; j++) {
//     //     integral += exp(lambda[j] * theta + log(gauss_weight[j] * omega_weight[j]));
//     //     printf("\tv1: exp(%f)\n", lambda[j] * theta + log(gauss_weight[j] * omega_weight[j]));
//     // }
//     //
//     // integral = log(integral);
//     //
//     // for(unsigned i=0; i<nint; i++) {
//     //     hull_estimate[i].observe(exp(lambda[i] * theta - integral));
//     //     printf("\tv2: exp(%f)\n", lambda[i] * theta - integral);
//     //     printf("value[%i] = %e\n", i, exp(lambda[i] * theta - integral));
//     // }
//
//
//     for(unsigned i=0; i<nint; i++) {
//         // holder for value
//         double value = 0.0;
//
//         vector<double> v(nint, 0.0);
//
//         // inner loop over all the points
//         for(unsigned j=0; j<nint; j++) {
//             // good for large distance between end points
//             value += exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
//         }
//
//         hull_estimate[i].observe(1.0 / value);
//     }
// }
//
// // update the weights
// void InfiniteSwitchSimulatedAppliedField::update_weights()
// {
//     // update the hull estimate
//     update_hull();
//
//     for(unsigned i=0; i<nint; i++) {
//         omega_weight[i] = (1.0 - time_step) * omega_weight[i]
//                             + time_step / hull_estimate[i].get_average();
//     }
//
//     // normalise the weights
//     normalize_weights();
// }
//
// void InfiniteSwitchSimulatedAppliedField::normalize_weights()
// {
//     double integral = 0.0;
//
//     for(unsigned i=0; i<nint; i++) {
//         integral += gauss_weight[i] * omega_weight[i];
//     }
//
//     for(auto& w : omega_weight) {
//         w /= integral;
//     }
// }
//
// // initialise the forces
// void InfiniteSwitchSimulatedAppliedField::initialize_force()
// {
//     // apply the force rescaling
//     apply_force_rescaling();
// }
//
// // get observable weights
// vector<double> InfiniteSwitchSimulatedAppliedField::get_observable_weights()
// {
//     // dereference values
//     double theta = get_collective();
//
//     vector<double> weights(nint, 1.0);
//
//     for(unsigned i=0; i<nint; i++) {
//
//         double integral = 0.0;
//
//         for(unsigned j=0; j<nint; j++) {
//             integral += exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
//         }
//
//         weights[i] /= hull_estimate[i].get_average() * integral;
//     }
//
//     return weights;
// }
