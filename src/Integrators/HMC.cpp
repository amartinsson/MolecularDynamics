#include "HMC.hpp"

using namespace::std;

// constructor
HMC::HMC(const double& time_step, const double& beta,
    System* system_pt, const int& seed)
        : Hamilton(system_pt), MetropolisHastings(beta, 0.0, system_pt, seed),
            time_step(time_step) {
                Hnp1 = 0.0;
                H = 0.0;
            }

// destructor
HMC::~HMC() {/* empty */}

// integrate
void HMC::integrate(Molecule* molecule_pt) {
    // save current positons
    log_current_position(molecule_pt);

    // solve for (q_np1, p_np1) using BAB and update hamiltoian
    Hnp1 = 0.0;
    hamiltoniain_forward(molecule_pt);

    // calculate change in energy
    double dE = Hnp1 - H;
    bool step = accept_reject(dE, proposal_ratio(molecule_pt));

    if(step) { // step was accepted
        // update hamiltoniain
        H = Hnp1;
    }
    else { // step was rejected
        for(auto& particle : molecule_pt->Particles) {
            // reset position
            particle.second->q = position[particle.second];
            // flip the momentum
            particle.second->p = momentum[particle.second].neg();
        }
        // reset the force
        Hamilton::System_pt->compute_force(molecule_pt);
    }
}

void HMC::log_current_position(Molecule* molecule_pt) {
    if(first_step) {
        // update the force
        Hamilton::System_pt->compute_force(molecule_pt);

        for(const auto& particle : molecule_pt->Particles) {
            // insert particle in position
            position.insert(make_pair(particle.second, particle.second->q));
            // update momentum
            propsal_momentum(*particle.second);
            momentum.insert(make_pair(particle.second, particle.second->p));
            // add to the hamiltoian
            H += calculate_hamiltonian(*particle.second);
        }

        // log that we have created lists
        first_step = false;
    }
    else {
        for(const auto& particle : molecule_pt->Particles) {
            // update position
            position[particle.second] = particle.second->q;
            // update momentum
            momentum[particle.second] = particle.second->p;
        }
    }
}

void HMC::propsal_momentum(Particle& particle) {
    // dereference beta
    double beta = MetropolisHastings::beta;

    // draw normal random number
    Vector N = Vector(particle.p.size(), MetropolisHastings::normal_gen);

    // draw new momentum
    particle.p = particle.m.inv() * N * beta;
}

double HMC::calculate_hamiltonian(const Particle& particle) {
    return 0.5 * particle.p.dot(particle.m.inv() * particle.p)
            + Hamilton::System_pt->compute_potential(particle.q);
}

void HMC::hamiltoniain_forward(Molecule* molecule_pt) {
    // loop over all the particles
    for(const auto& particle : molecule_pt->Particles) {
        // momentum
        Hamilton::B(*particle.second, 0.5 * time_step);
        // position
        Hamilton::A(*particle.second, time_step);
    }

    // solve for the force
    MetropolisHastings::system->compute_force(molecule_pt);

    // Do final B step
    for(const auto& particle : molecule_pt->Particles) {
        // momentum
        Hamilton::B(*particle.second, 0.5 * time_step);
        // add to the new Hamiltoniain
        Hnp1 += calculate_hamiltonian(*particle.second);
    }
}
