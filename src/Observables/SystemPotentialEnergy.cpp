#include "SystemPotentialEnergy.hpp"

// constructor
SystemPotentialEnergy::SystemPotentialEnergy(const Molecule* molecule_pt,
    const int& recf, const int& rect)
        : SystemEnergy(molecule_pt, recf, rect) {/* empty */ }

// update energy
void SystemPotentialEnergy::update() {
    if(SystemEnergy::recStep()) {
        // get potential energy and Kinetic energy
        double V = SystemEnergy::molecule_pt->potential();

        // add this to observation
        SystemEnergy::Energy->observe(V);
    }
}
