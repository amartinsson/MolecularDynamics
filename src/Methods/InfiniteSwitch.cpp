#include "InfiniteSwitch.hpp"

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
// constructor
InfiniteSwitch::InfiniteSwitch(Molecule* molecule_pt,
    const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step)
         : molecule_pt(molecule_pt), potential_pt(&molecule_pt->potential()),
            lambda_min(lambda_min), lambda_max(lambda_max), nint(nint),
                time_step(time_step), beta(molecule_pt->beta())
{
    // initialize vectors
    gauss_weight = new double[nint];

    lambda = new double[nint];

    hull_estimate = new AverageObservable[nint];

    omega_weight = std::vector<double>(nint, 1.0);

    // initialize integral calculation
    legendre_compute_glr(nint, lambda, gauss_weight);

    rescale(lambda_min, lambda_max, nint, lambda, gauss_weight);
}

// destructor
InfiniteSwitch::~InfiniteSwitch()
{
    for(unsigned i=0; i<nint; i++) {
        delete &gauss_weight[i];
        delete &lambda[i];
        delete &hull_estimate[i];
    }
}

// get the thermal rescaling factor
void InfiniteSwitch::apply_force_rescaling()
{
    // rescaling
    double lambda_bar = calculate_lambda_bar();

    double kT = 1.0 / beta;

    for(auto& particle : molecule_pt->Particles) {
        // get the gradient of the collective variable
        Vector grad = *grad_collective_pt[particle.second];

        // add to the rescaling of the force
        particle.second->f -= grad * kT * lambda_bar;

        if(particle.second->rigid_body() == true) {
            printf("ERROR: Infinite Switch with rigid_body not implimented\n");
            exit(-1);
        }
    }
}

// must have private calculate force rescaling
double InfiniteSwitch::calculate_lambda_bar()
{
    // update the weights at current position
    update_weights();

    // do integrals
    double lambda_bar_lower = 0.0;
    double lambda_bar_upper = 0.0;

    double V = *potential_pt;
    double theta = *collective_pt;

    for(unsigned i=0; i<nint; i++) {
        // double measure = omega_weight[i] * exp(-beta * V + lambda[i] * theta);
        double measure = omega_weight[i] * exp(lambda[i] * theta);

        lambda_bar_upper += gauss_weight[i] * lambda[i] * measure;
        lambda_bar_lower += gauss_weight[i] * measure;
    }

    return lambda_bar_upper / lambda_bar_lower;
}

// update the hull estimate
void InfiniteSwitch::update_hull()
{
    // dereference values
    double theta = *collective_pt;

    // do the integral first
    double integral = 0.0;

    for(unsigned i=0; i<nint; i++) {
        integral += gauss_weight[i] * omega_weight[i]
                        *  exp(lambda[i] * theta);
    }

    // update the hull estimate
    for(unsigned i=0; i<nint; i++) {
        double value = exp(lambda[i] * theta) / integral;

        hull_estimate[i].observe(value);
    }
}

// update the weights
void InfiniteSwitch::update_weights()
{
    // update the hull estimate
    update_hull();

    for(unsigned i=0; i<nint; i++) {
        omega_weight[i] = (1.0 - time_step) * omega_weight[i]
                            + time_step / hull_estimate[i].get_average();
    }

    // normalise the weights
    normalize_weights();
}

void InfiniteSwitch::normalize_weights()
{
    double integral = 0.0;

    for(unsigned i=0; i<nint; i++) {
        integral += gauss_weight[i] * omega_weight[i];
    }

    for(auto& w : omega_weight) {
        w /= integral;
    }
}

// initialise the forces
void InfiniteSwitch::initialize_force()
{
    // apply the force rescaling
    apply_force_rescaling();
}

// get observable weights
vector<double> InfiniteSwitch::get_observable_weights()
{
    vector<double> weights(nint, 0.0);

    for(unsigned i=0; i<nint; i++) {
        weights[i] = hull_estimate[i].get_average();
    }

    return weights;
}
