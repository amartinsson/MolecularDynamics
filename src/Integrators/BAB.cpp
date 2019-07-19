#include "BAB.hpp"

// constructor
BAB::BAB(const double& time_step, System* system_pt)
    : Hamilton(system_pt), time_step(time_step) {/* empty */}

// integrator
void BAB::integrate(Molecule* molecule_pt) {
    // loop over all the particles
    for(const auto& particle : molecule_pt->Particles) {
        // momentum
        Hamilton::B(*particle.second, 0.5 * time_step);
        // position
        Hamilton::A(*particle.second, time_step);
    }

    // solve for the force
    Hamilton::System_pt->compute_force(molecule_pt);

    // Do final B step
    for(const auto& particle : molecule_pt->Particles) {
        // momentum
        Hamilton::B(*particle.second, 0.5 * time_step);
    }
}
