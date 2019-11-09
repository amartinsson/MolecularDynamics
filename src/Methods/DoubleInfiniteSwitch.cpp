#include "DoubleInfiniteSwitch.hpp"

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
// constructor
DoubleInfiniteSwitch::DoubleInfiniteSwitch(Molecule* molecule_pt,
    const double& lambda_one_min, const double& lambda_one_max,
        const double& lambda_two_min, const double& lambda_two_max,
            const unsigned& nOneint, const unsigned& nTwoint,
                const double& time_step)
    : molecule_pt(molecule_pt), potential_pt(&molecule_pt->potential()),
        lambda_one_min(lambda_one_min), lambda_one_max(lambda_one_max),
            lambda_two_min(lambda_two_min), lambda_two_max(lambda_two_max),
                nOneint(nOneint), nTwoint(nTwoint), time_step(time_step),
                    beta(molecule_pt->beta()), omega_weight(nOneint, nTwoint)
{
    // initialize vectors
    gauss_one_weight = new double[nOneint];
    gauss_two_weight = new double[nTwoint];

    lambda_one = new double[nOneint];
    lambda_two = new double[nTwoint];

    for(unsigned i=0; i<nOneint; i++) {
        vector<AverageObservable*> obs;

        for(unsigned j=0; j<nTwoint; j++) {
            obs.push_back(new AverageObservable());

            // set the omega weights to one
            omega_weight(i, j) = 1.0;
        }

        hull_estimate.insert(make_pair(i, obs));
    }

    // initialize integral calculation
    legendre_compute_glr(nOneint, lambda_one, gauss_one_weight);
    legendre_compute_glr(nTwoint, lambda_two, gauss_two_weight);

    rescale(lambda_one_min, lambda_one_max, nOneint, lambda_one,
        gauss_one_weight);
    rescale(lambda_two_min, lambda_two_max, nTwoint, lambda_two,
        gauss_two_weight);
}

// destructor
DoubleInfiniteSwitch::~DoubleInfiniteSwitch()
{
    // for(unsigned i=0; i<nint; i++) {
    //     delete &gauss_weight[i];
    //     delete &lambda[i];
    //     delete &hull_estimate[i];
    // }
}

// get the thermal rescaling factor
void DoubleInfiniteSwitch::apply_force_rescaling()
{
    // rescaling
    Vector lambda_bar = calculate_lambda_bar();

    double kT = 1.0 / beta;

    for(auto& particle : molecule_pt->Particles) {
        // get the gradient of the collective variable
        // Vector grad_one = *grad_collective_one_pt[particle.second];
        // Vector grad_two = *grad_collective_two_pt[particle.second];
        // get the gradient of the collective variable
        Vector grad_one = get_collective_grad_one(particle.second);
        Vector grad_two = get_collective_grad_two(particle.second);

        // // add to the rescaling of the force
        // particle.second->f -= grad_one * kT * lambda_bar(0)
        //                         + grad_two * kT * lambda_bar(1);
        particle.second->f += grad_one * kT * lambda_bar(0)
                                + grad_two * kT * lambda_bar(1);

        if(particle.second->rigid_body() == true) {
            printf("ERROR: Double Infinite Switch with rigid_body not implimented\n");
            exit(-1);
        }
    }
}

// must have private calculate force rescaling
Vector DoubleInfiniteSwitch::calculate_lambda_bar()
{
    // update the weights at current position
    update_weights();

    // do integrals
    double lambda_bar_lower = 0.0;

    double lambda_bar_one_upper = 0.0;
    double lambda_bar_two_upper = 0.0;

    double V = *potential_pt;

    // double theta_one = *collective_one_pt;
    // double theta_two = *collective_two_pt;
    double theta_one = get_collective_one();
    double theta_two = get_collective_two();

#pragma omp parallel for collapse(2) schedule(static) firstprivate(beta, V, theta_one, theta_two, nOneint, nTwoint)\
    reduction(+:lambda_bar_lower,lambda_bar_one_upper, lambda_bar_two_upper)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {
            // either of these will work, first should increase stability
            // of exp calculation
            double measure = omega_weight(i, j) * exp(-beta * V + lambda_one[i] * theta_one + lambda_two[j] * theta_two);

            lambda_bar_lower += gauss_one_weight[i] * gauss_two_weight[j]
                                        * measure;

            lambda_bar_one_upper += gauss_one_weight[i] * gauss_two_weight[j]
                                        * lambda_one[i] * measure;

            lambda_bar_two_upper += gauss_one_weight[i] * gauss_two_weight[j]
                                        * lambda_two[j] * measure;
        }
    }

    // vector to return
    Vector lambda_bar(2);

    lambda_bar(0) = lambda_bar_one_upper / lambda_bar_lower;
    lambda_bar(1) = lambda_bar_two_upper / lambda_bar_lower;

    return lambda_bar;
}

// update the hull estimate
void DoubleInfiniteSwitch::update_hull()
{
    // dereference values
    // double theta_one = *collective_one_pt;
    // double theta_two = *collective_two_pt;
    double theta_one = get_collective_one();
    double theta_two = get_collective_two();

#pragma omp parallel for collapse(2) schedule(static) firstprivate(theta_one, theta_two, nOneint, nTwoint)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {

            // holder for value
            double value = 0.0;

            // inner loop over all the points
            for(unsigned k=0; k<nOneint; k++) {
                for(unsigned kk=0; kk<nTwoint; kk++) {

                    value += gauss_one_weight[k] * gauss_two_weight[kk] * omega_weight(k, kk)
                    * exp((lambda_one[k] - lambda_one[i]) * theta_one + (lambda_two[kk] - lambda_two[j]) * theta_two);
                }
            }

            hull_estimate[i][j]->observe(1.0 / value);
        }
    }
}

// update the weights
void DoubleInfiniteSwitch::update_weights()
{
    // update the hull estimate
    update_hull();

#pragma omp parallel for collapse(2) schedule(static) firstprivate(time_step, nOneint, nTwoint)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {
            omega_weight(i,j) = (1.0 - time_step) * omega_weight(i, j)
                + time_step / hull_estimate[i][j]->get_average();
        }
    }

    // normalise the weights
    normalize_weights();
}

void DoubleInfiniteSwitch::normalize_weights()
{
    double integral = 0.0;

#pragma omp parallel for collapse(2) schedule(static) reduction(+:integral) firstprivate(nOneint, nTwoint)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {
            integral += gauss_one_weight[i] * gauss_two_weight[j] * omega_weight(i,j);
        }
    }
#pragma omp parallel for collapse(2) schedule(static) firstprivate(nOneint, nTwoint)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {
            omega_weight(i,j) /= integral;
        }
    }

}

// initialise the forces
void DoubleInfiniteSwitch::initialize_force()
{
    // apply the force rescaling
    apply_force_rescaling();
}

// get observable weights
Matrix DoubleInfiniteSwitch::get_observable_weights()
{
    // dereference values
    // double theta_one = *collective_one_pt;
    // double theta_two = *collective_two_pt;
    double theta_one = get_collective_one();
    double theta_two = get_collective_two();

    Matrix weights(nOneint, nTwoint);

#pragma omp parallel for collapse(2) schedule(static) firstprivate(theta_one, theta_two, nOneint, nTwoint)
    for(unsigned i=0; i<nOneint; i++) {
        for(unsigned j=0; j<nTwoint; j++) {

            double integral = 0.0;

            // inner loop over all the points
            for(unsigned k=0; k<nOneint; k++) {
                for(unsigned kk=0; kk<nTwoint; kk++) {

                    integral += gauss_one_weight[k] * gauss_two_weight[kk] * omega_weight(k, kk)
                    * exp((lambda_one[k] - lambda_one[i]) * theta_one + (lambda_two[kk] - lambda_two[j]) * theta_two);
                }
            }

                weights(i,j) = 1.0 / (hull_estimate[i][j]->get_average() * integral);
        }
    }

    return weights;
}


double DoubleInfiniteSwitch::get_mid_lambda_one() {
    return lambda_one[nOneint / 2];
}

double DoubleInfiniteSwitch::get_mid_lambda_two() {
    return lambda_two[nTwoint / 2];
}
