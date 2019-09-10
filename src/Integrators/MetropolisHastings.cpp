#include "MetropolisHastings.hpp"

using namespace::std;

MetropolisHastings::MetropolisHastings(const double& beta, const double& sigma,
    System* system, const int& seed)
        : normal_gen(0.0, 1.0, seed), beta(beta), sigma(sigma), system(system),
            uniform_gen(0.0, 1.0, seed+1), first_step(true) {
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
    // for first step make sure map exists
    if(first_step)
        create_position_map(molecule_pt);

    // loop over all the particles
    Vnp1 = 0.0;
    for(auto& particle : molecule_pt->Particles) {
        // generate new proposal
        qprime[particle.second] = proposal_move(particle.second->q);

        // calclate the potential at this step
        Vnp1 += system->compute_potential(qprime[particle.second]);
    }

    // calculate the difference in energy
    double dE = Vnp1 - V;
    bool step = accept_reject(dE, proposal_ratio(molecule_pt));

    // take action based on outcome of step
    if(step) { // step accepted
        // update the potential
        V = Vnp1;
        // update the positions
        for(auto& particle : molecule_pt->Particles) {
            particle.second->q = qprime[particle.second];
        }
        // cumpute the dummy variables
        system->compute_force(molecule_pt);
    }
    else { // step rejected
        // do nothing as nothing changes
    }
}

Vector MetropolisHastings::proposal_move(const Vector& q) {
    // make holder for new vector
    Vector qprime = Vector(q.size());
    Vector N = Vector(q.size(), normal_gen);

    // generate new position
    qprime = q + N * sigma;

    // return the newly suggested position
    return qprime;
}


double MetropolisHastings::proposal_ratio(Molecule* molecult_pt) {
    // in this case it's symmetric so let's not worry
    return 1.0;
}

void MetropolisHastings::create_position_map(Molecule* molecule_pt) {
    // update the force for dummy variables
    system->compute_force(molecule_pt);

    for(const auto& particle : molecule_pt->Particles) {
        // insert particle in position
        qprime.insert(make_pair(particle.second, particle.second->q));

        // add to the potential
        V += system->compute_potential(qprime[particle.second]);
    }

    // log that we have created lists
    first_step = false;
}

// make an accept reject evaulation with acceptance a
bool MetropolisHastings::accept_reject(const double& dE,
    const double& proposal_ratio) {
        // return boolean
        bool accepted = false;

        // calculate acceptance
        double a = exp(-beta * dE) * proposal_ratio;

        // make a uniform
        double u = uniform_gen();

        if(u <= a) {
            // log the step as accepted
            acceptance->observe(1.0);
            // change boolean to accepted
            accepted = true;
        }
        else {
            // log the step as rejected
            acceptance->observe(0.0);
        }

        // return boolean if it was accepted
        return accepted;
    }
