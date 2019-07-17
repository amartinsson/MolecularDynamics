#include "MetropolisHastings.hpp"

using namespace::std;

MetropolisHastings::MetropolisHastings(const double& beta, const double& sigma,
    System* system, const int& seed)
        : normal_gen(0.0, 1.0, seed), beta(beta), sigma(sigma), system(system),
            uniform_gen(0.0, 1.0, seed+1) {
    // make a new average object to keep track of the acceptance ratio
    acceptance = new AverageObservable();
}

MetropolisHastings::~MetropolisHastings() {

    // log the ratio to prompt
    printf("\nAcceptance Ratio: %2.2f\n", acceptance->get_average() * 100.0);

    // delete the object
    delete acceptance;
}

// integrator
void MetropolisHastings::integrate(Molecule* molecule_pt) {
    // reset the potential
    molecule_pt->potential() = 0.0;

    // loop over all the particles
    for(auto& particle : molecule_pt->Particles) {
        // generate new proposal
        Vector qprime = proposal_move(particle.second->q);

        // calculate acceptance ratio
        double a = acceptance_ratio(particle.second->q, qprime);

        // make a uniform
        double u = uniform_gen();

        if(u <= a) {
            // accept the step
            particle.second->q = qprime;
            // add to the potential
            molecule_pt->potential() += system->compute_potential(particle.second->q);
            // log the step as accepted
            acceptance->observe(1.0);
        }
        else {
            // log the step as rejected
            acceptance->observe(0.0);
            // add to the potential
            molecule_pt->potential() += system->compute_potential(particle.second->q);
        }
    }
}

Vector MetropolisHastings::proposal_move(const Vector& q) {
    // make holder for new vector
    Vector qprime = Vector(q.size());
    Vector N = Vector(q.size(), normal_gen);

    // generate new position
    qprime = q + N * pow(sigma, 2.0);

    // return the newly suggested position
    return qprime;
}


double MetropolisHastings::proposal_ratio(const Vector& q,
                                          const Vector& qprime) {
    // in this case it's symmetric so let's not worry
    return 1.0;
}

double MetropolisHastings::acceptance_ratio(const Vector& q,
    const Vector& qprime) {
        // delta energy
        double dE = system->compute_potential(qprime)
                    - system->compute_potential(q);

        // holder for acceptance
        return exp(-beta * dE) * proposal_ratio(q, qprime);
}
