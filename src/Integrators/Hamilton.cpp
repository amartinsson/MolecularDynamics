#include "Hamilton.hpp"

// constructor
Hamilton::Hamilton(System* system_pt)
    : System_pt(system_pt) {/* empty */}

Hamilton::~Hamilton() {/* empty */}

// A step
void Hamilton::A(Particle& particle, const double& h) {
    // update position
    particle.q += particle.m.inv() * particle.p * h;
}

// B step
void Hamilton::B(Particle& particle, const double& h) {
    // update momentum
    particle.p += particle.f * h;
}

// // A step
// void Hamilton::A(Vector& q, const Vector& m, const Vector& p, const double& h) {
//     // update position
//     q += m.inv() * p * h;
// }
//
// // B step
// void Hamilton::B(Vector& p, const Vector& f, const double& h) {
//     // update momentum
//     p += f * h;
// }

// integrator
void Hamilton::integrate(Molecule* molecule_pt) {
    std::cout << "Error: Hamilton class must implement integrate() function\n";
    exit(-1);
}
