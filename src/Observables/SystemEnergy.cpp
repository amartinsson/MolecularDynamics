#include "SystemEnergy.hpp"

using namespace::std;

// constructor
SystemEnergy::SystemEnergy(const Molecule* molecule_pt, const int& recf,
                           const int& rect)
            : SystemObservable(recf, rect), molecule_pt(molecule_pt) {
                // make new energy object
                Energy = new AverageObservable();
            }

// destructor
SystemEnergy::~SystemEnergy() {
    delete Energy;
}

// instant
double SystemEnergy::get_instant() {
    return Energy->get_instant();
}

// average
double SystemEnergy::get_average() {
    return Energy->get_average();
}

// update function
void SystemEnergy::update() {
    if(recStep()) {
        // get potential energy and Kinetic energy
        double V = molecule_pt->potential();
        double K = this->kinetic();

        // add this to observation
        Energy->observe(K + V);
    }
}

// calculate the kinetic energy
double SystemEnergy::kinetic() {
    // holder for kinetic energy
    double K = 0.0;

    // loop over the number of particles
    for(auto& particle : molecule_pt->Particles) {
        K += particle.second->p.dot(particle.second->m.inv() * particle.second->p);
    }

    // divide by 2 and return
    return 0.5 * K;
}
