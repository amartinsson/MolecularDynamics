#include "AdaptiveBiasingForce.hpp"
#include <algorithm>
#include <numeric>

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
// constructor
AdaptiveBiasingForce::AdaptiveBiasingForce(Molecule* molecule_pt,
    const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
            const unsigned& threshold)
         : molecule_pt(molecule_pt), potential_pt(&molecule_pt->potential()),
            lambda_min(lambda_min), lambda_max(lambda_max), nint(nint),
                time_step(time_step), beta(molecule_pt->beta()),
                    StepCount(0), Threshold(threshold), ilocator(HistObservable(lambda_min, lambda_max, nint))
{
    partition_estimate = new AverageObservable[nint];

    force_estimate = new AverageObservable[nint];
}

// destructor
AdaptiveBiasingForce::~AdaptiveBiasingForce()
{
    // print the pdf of the biasing force
    for(unsigned i=0; i<nint; i++) {
        double x = ilocator.get_lower_bin(i);
        double abf = force_estimate[i].get_average();
        int n = force_estimate[i].get_nobservations();
        printf("%1.3f, %1.3f, %.0d\n", x, abf, n);
    }

    // for(unsigned i=0; i<nint; i++) {
    //     delete &partition_estimate[i];
    //     delete &force_estimate[i];
    // }
}

// get the thermal rescaling factor
void AdaptiveBiasingForce::apply_force_rescaling()
{
    // update the weights at current position
    // update_weights();

    // get the index for which we are evaluating
    int index = ilocator.find_index(get_collective());

    // get the number of observations at this level
    int nobs = force_estimate[index].get_nobservations();

    if(nobs > Threshold) {
        // if we have recorded enough points at this level
        // apply the force rescaling
        update_adaptive_force(index);

        // apply the rescaling
        apply_force_rescaling(index);
    }
    else {
        // update the force for the correct index
        update_adaptive_force(index);
    }
}

void AdaptiveBiasingForce::update_adaptive_force(const int& index)
{
    // holder for the force
    double abf = 0.0;
    double V = *potential_pt;

    // get the particle
    Particle particle = molecule_pt->particle(0);

    // get the original force on the particle
    Vector grad = particle.f.neg();

    // get the gradient of the collective varialbe
    Vector grad_xi = get_collective_grad(&particle);

    // get the laplacian of the collective varialbe
    double lapl_xi = get_collective_laplace(&particle);

    // dot product of the gradient
    double grad_grad = grad_xi.dot(grad_xi);

    // calculate the adaptive biasing force
    abf = -grad.dot(grad_xi) / grad_grad - 1.0 / beta * lapl_xi / grad_grad;

    // printf("x = %f, V = %f, f = %f, dxi = %f, |grad|^2 = %f, ddxi = %f, abf = %f\n",
    // particle.q(0), V, particle.f(0), grad_xi(0), grad_grad, lapl_xi, abf);

    // observe the force
    force_estimate[index].observe(abf);

    // observe the partition estimate
    partition_estimate[index].observe(exp(-beta * V));
}

// get the thermal rescaling factor
void AdaptiveBiasingForce::apply_force_rescaling(const int& index)
{
    // printf("lambda_bar = %f\n", lambda_bar);
    double kT = 1.0 / beta;

    // get the particle
    Particle* particle = &molecule_pt->particle(0);

    // get the gradient of the collective variable
    Vector grad = get_collective_grad(particle);

    double Gprime = -force_estimate[index].get_average();

    // printf("original force %f ", particle.f(0));
    // add to the rescaling of the force
    particle->f += grad * kT * Gprime;

    // printf("weighted force %f\n", particle.f(0));

}

// // must have private calculate force rescaling
// double AdaptiveBiasingForce::calculate_lambda_bar()
// {
//     // do integrals
//     double lambda_bar_lower = 0.0;
//     double lambda_bar_upper = 0.0;
//
//     // double V = *potential_pt;
//     //
//     // // add to the potential if we're using npt
//     // if(with_npt) {
//     //     // V += *pressure_pt * (*cell_pt).det();
//     //     V += npt_pt->get_target_pressure() * npt_pt->S.det();
//     // }
//
//     double theta = get_collective();
//
//     vector<double> lbar_u(nint, 0.0);
//     vector<double> lbar_l(nint, 0.0);
//
//     for(unsigned i=0; i<nint; i++) {
//         // version one
//         // // either of these will work, first should increase stability
//         // // of exp calculation
//         // // double measure = omega_weight[i] * exp(-beta * V + lambda[i] * theta);
//         // double measure = omega_weight[i] * exp(lambda[i] * theta);
//         //
//         //
//         // // // // printf("\t exp(%f)\n", -beta * V + lambda[i] * theta);
//         // // printf("\t exp(%f) and log(%f) vs exp(%f)\n",
//         // // lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])),
//         // // omega_weight[i] * gauss_weight[i] * lambda[i],
//         // // lambda[i] * theta);
//         // //
//         // printf("%i exp[%f]\n", i, lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
//         // lambda_bar_upper += gauss_weight[i] * lambda[i] * measure;
//         // lambda_bar_lower += gauss_weight[i] * measure;
//
//         // lambda_bar_upper += exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
//         // lambda_bar_lower += exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
//
//         lbar_u[i] = exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
//         lbar_l[i] = exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
//     }
//
//     // sort in ascending order
//     sort(lbar_u.begin(), lbar_u.end());
//     sort(lbar_l.begin(), lbar_l.end());
//
//     double su = 0.0;
//     double sl = 0.0;
//
//     accumulate(lbar_u.begin(), lbar_u.end(), su);
//     accumulate(lbar_l.begin(), lbar_l.end(), sl);
//
//     double upper = log(lbar_u[nint-1]) + log(1.0 + su / lbar_u[nint-1]);
//     double lower = log(lbar_l[nint-1]) + log(1.0 + sl / lbar_l[nint-1]);
//
//
//     // exit(-1);
//
//     // printf("rescale: %e, (trick) %e\n", lambda_bar_upper / lambda_bar_lower, exp(upper - lower));
//     // printf("(trick) %e\n",exp(upper - lower));
//
//     // return lambda_bar_upper / lambda_bar_lower;
//     return exp(upper - lower);
// }
//
// // update the hull estimate
// void AdaptiveBiasingForce::update_hull()
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
//         // double value = 0.0;
//
//         vector<double> v(nint, 0.0);
//
//         // inner loop over all the points
//         for(unsigned j=0; j<nint; j++) {
//
//
//             // good for large distance between end points
//             // value += exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
//
//             // value += gauss_weight[j] * omega_weight[j] * exp((lambda[j] - lambda[i]) * theta);
//             // printf("\t %i exp[%f]\n", i, (lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
//
//             v[j] = exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
//         }
//
//         // sort in ascending order
//         sort(v.begin(), v.end());
//
//         double s = 0.0;
//
//         accumulate(v.begin(), v.end(), s);
//
//         double value2 = log(v[nint-1]) + log(1.0 + s/v[nint-1]);
//
//         // printf("value = %e, value(trick) = %e\n", value, exp(value2));
//         // printf("value(trick) = %e\n", exp(-value2));
//         // printf("---- value[%d] %e\n", i, 1.0 / value);
//         // hull_estimate[i].observe(1.0 / value);
//         hull_estimate[i].observe(exp(-value2));
//     }
//     // exit(-1);
// }
//
// // update the weights
// void AdaptiveBiasingForce::update_weights()
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
// void AdaptiveBiasingForce::normalize_weights()
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
// void AdaptiveBiasingForce::initialize_force()
// {
//     // apply the force rescaling
//     apply_force_rescaling();
// }
//
// // initialize with npt
// void AdaptiveBiasingForce::initialize_with_npt(NptGrid* npt_grid_pt)
// {
//     // set to use npt
//     with_npt = true;
//
//     // set the pointer to the npt grid
//     npt_pt = npt_grid_pt;
// }
//
// // get observable weights
// vector<double> AdaptiveBiasingForce::get_observable_weights()
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
