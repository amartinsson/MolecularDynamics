#include "SystemVirialTemperature.hpp"

// constructor
SystemVirialTemperature::SystemVirialTemperature(
    Molecule* molecule_pt, const int& recf, const int& rect)
        : SystemTemperature(molecule_pt, recf, rect) {}

// destructor
SystemVirialTemperature::~SystemVirialTemperature() {}

// update function
void SystemVirialTemperature::update()
{
    if(recStep()) {
        // update the temperature
        update_temperature();
    }
}

// calculate the Virial temperature
void SystemVirialTemperature::update_temperature() {
    // make temporary temperature
    double T = 0.0;

    double dim = (double) system->dim();
    double N = (double) system->nparticle();

    for(const auto& particle : system->Particles) {
        T = -1.0 / (N * dim) * particle.second->q.dot(particle.second->f)
            + (N - 1.0) / N * T;
    }

    // add to nabla square
    Temp->observe(T);
}
