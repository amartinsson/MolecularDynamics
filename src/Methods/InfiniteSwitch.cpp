#include "InfiniteSwitch.hpp"
#include <algorithm>
#include <numeric>

using namespace::std;

/******************************************************************************
                           Infinite Switch Simulation
 *****************************************************************************/
// constructor
InfiniteSwitch::InfiniteSwitch(Molecule* molecule_pt,
    const double& lambda_min, const double& lambda_max,
        const unsigned& nint, const double& time_step,
            const unsigned& threshold)
         : molecule_pt(molecule_pt), potential_pt(&molecule_pt->potential()),
            lambda_min(lambda_min), lambda_max(lambda_max), nint(nint),
                time_step(time_step), beta(molecule_pt->beta()),
                    StepCount(0), Threshold(threshold)
{
    // initialize vectors
    gauss_weight = new double[nint];

    lambda = new double[nint];

    hull_estimate = new AverageObservable[nint];

    omega_weight = std::vector<double>(nint, 1.0);

    // initialize integral calculation
    legendre_compute_glr(nint, lambda, gauss_weight);

    rescale(lambda_min, lambda_max, nint, lambda, gauss_weight);

    for(unsigned i=0; i<nint; i++)
        printf("l[%i] = %f\n", i, lambda[i]);
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
    // update the weights at current position
    update_weights();

    if(StepCount > Threshold) {
        // apply the rescaling
        apply_force_rescaling_all();
    }
    else {
        StepCount++;

    }
}

// get the thermal rescaling factor
void InfiniteSwitch::apply_force_rescaling_all()
{
    // rescaling
    double lambda_bar = calculate_lambda_bar();

    // printf("lambda_bar = %f\n", lambda_bar);
    double kT = 1.0 / beta;

    for(auto& particle : molecule_pt->Particles) {
        // get the gradient of the collective variable
        Vector grad = get_collective_grad(particle.second);

        // add to the rescaling of the force
        particle.second->f += grad * kT * lambda_bar;

        if(particle.second->rigid_body() == true) {
            particle.second->tau += get_collective_grad_rot(particle.second) *
                                        kT * lambda_bar;
        }
    }

    if(with_npt) {
        // get the dimension
        unsigned dim = molecule_pt->dim();

        // get collective virial grad
        Matrix grad = get_collective_virial_grad();

        npt_pt->virial += grad * kT * lambda_bar;

        // // rescale the virial function
        // for(unsigned i=0; i<dim; i++) {
        //     for(unsigned j=i; j<dim; j++) {
        //         npt_pt->virial(i, j) += grad(i, j) * kT * lambda_bar;
        //     }
        // }
    }
}

// must have private calculate force rescaling
double InfiniteSwitch::calculate_lambda_bar()
{

    // double V = *potential_pt;
    //
    // // add to the potential if we're using npt
    // if(with_npt) {
    //     // V += *pressure_pt * (*cell_pt).det();
    //     V += npt_pt->get_target_pressure() * npt_pt->S.det();
    // }

    double theta = get_collective();
    double lambda_bar = 0.0;


    if(lambda[0] > 0.0) {

        vector<double> lbar_u(nint, 0.0);
        vector<double> lbar_l(nint, 0.0);

        for(unsigned i=0; i<nint; i++) {
            // version one
            // // either of these will work, first should increase stability
            // // of exp calculation
            // // double measure = omega_weight[i] * exp(-beta * V + lambda[i] * theta);
            // double measure = omega_weight[i] * exp(lambda[i] * theta);
            //
            //
            // // // // printf("\t exp(%f)\n", -beta * V + lambda[i] * theta);
            // // printf("\t exp(%f) and log(%f) vs exp(%f)\n",
            // // lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])),
            // // omega_weight[i] * gauss_weight[i] * lambda[i],
            // // lambda[i] * theta);
            // //
            // printf("%i exp[%f]\n", i, lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
            // lambda_bar_upper += gauss_weight[i] * lambda[i] * measure;
            // lambda_bar_lower += gauss_weight[i] * measure;

            // lambda_bar_upper += exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
            // lambda_bar_lower += exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));

            lbar_u[i] = exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i] * lambda[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));
            lbar_l[i] = exp(lambda[i] * theta + log(omega_weight[i] * gauss_weight[i]) - (lambda[0] * theta + log(omega_weight[0] * gauss_weight[0] * lambda[0])));

            // printf("lambda %e theta %e omega_weight %e gauss_weight %e\n", lambda[i], theta, omega_weight[i], gauss_weight[i]);
            // printf("log %e\n", omega_weight[i] * gauss_weight[i] * lambda[i]);
            // printf("lbar_u[%i] = %e\n", i, lbar_u[i]);
            // printf("lbar_l[%i] = %e\n", i, lbar_l[i]);
        }

        // sort in ascending order
        sort(lbar_u.begin(), lbar_u.end());
        sort(lbar_l.begin(), lbar_l.end());

        double su = 0.0;
        double sl = 0.0;

        accumulate(lbar_u.begin(), lbar_u.end(), su);
        accumulate(lbar_l.begin(), lbar_l.end(), sl);

        double upper = log(lbar_u[nint-1]) + log(1.0 + su / lbar_u[nint-1]);
        double lower = log(lbar_l[nint-1]) + log(1.0 + sl / lbar_l[nint-1]);

        lambda_bar = exp(upper - lower);
    }
    else {
        // do integrals
        double lambda_bar_lower = 0.0;
        double lambda_bar_upper = 0.0;

        for(unsigned i=0; i<nint; i++) {
            // calculate the measure
            double measure = omega_weight[i] * exp(lambda[i] * theta);

            // calculate upper and lower values
            lambda_bar_upper += gauss_weight[i] * lambda[i] * measure;
            lambda_bar_lower += gauss_weight[i] * measure;
        }

        lambda_bar = lambda_bar_upper / lambda_bar_lower;
    }



    // exit(-1);

    // printf("rescale: %e, (trick) %e\n", lambda_bar_upper / lambda_bar_lower, exp(upper - lower));
    // printf("(trick) %e\n",exp(upper - lower));

    // printf("returning lambda bar %e\n", lambda_bar);

    // return lambda_bar_upper / lambda_bar_lower;
    return lambda_bar;
}

// update the hull estimate
void InfiniteSwitch::update_hull()
{
    // dereference values
    double theta = get_collective();

    // double integral = 0.0;
    //
    // for(unsigned j=0; j<nint; j++) {
    //     integral += exp(lambda[j] * theta + log(gauss_weight[j] * omega_weight[j]));
    //     printf("\tv1: exp(%f)\n", lambda[j] * theta + log(gauss_weight[j] * omega_weight[j]));
    // }
    //
    // integral = log(integral);
    //
    // for(unsigned i=0; i<nint; i++) {
    //     hull_estimate[i].observe(exp(lambda[i] * theta - integral));
    //     printf("\tv2: exp(%f)\n", lambda[i] * theta - integral);
    //     printf("value[%i] = %e\n", i, exp(lambda[i] * theta - integral));
    // }


    for(unsigned i=0; i<nint; i++) {
        // holder for value
        // double value = 0.0;

        vector<double> v(nint, 0.0);

        // inner loop over all the points
        for(unsigned j=0; j<nint; j++) {


            // good for large distance between end points
            // value += exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));

            // value += gauss_weight[j] * omega_weight[j] * exp((lambda[j] - lambda[i]) * theta);
            // printf("\t %i exp[%f]\n", i, (lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));

            v[j] = exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
            // printf("v[%d] = %e\n", j, v[j]);
        }

        // sort in ascending order
        sort(v.begin(), v.end());

        double s = 0.0;

        accumulate(v.begin(), v.end(), s);

        double value2 = log(v[nint-1]) + log(1.0 + s/v[nint-1]);

        // printf("value = %e, value(trick) = %e\n", value, exp(value2));
        // printf("value(trick) = %e\n", exp(-value2));
        // printf("---- value[%d] %e\n", i, 1.0 / value);
        // hull_estimate[i].observe(1.0 / value);
        hull_estimate[i].observe(exp(-value2));
    }
    // exit(-1);
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

// initialize with npt
void InfiniteSwitch::initialize_with_npt(NptGrid* npt_grid_pt)
{
    // set to use npt
    with_npt = true;

    // set the pointer to the npt grid
    npt_pt = npt_grid_pt;
}

// get observable weights
vector<double> InfiniteSwitch::get_observable_weights()
{
    // dereference values
    double theta = get_collective();

    vector<double> weights(nint, 1.0);

    for(unsigned i=0; i<nint; i++) {

        double integral = 0.0;

        for(unsigned j=0; j<nint; j++) {
            integral += exp((lambda[j] - lambda[i]) * theta + log(gauss_weight[j] * omega_weight[j]));
        }

        weights[i] /= hull_estimate[i].get_average() * integral;
    }

    return weights;
}
