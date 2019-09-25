#include "DoubleWell.hpp"

DoubleWell::DoubleWell(const double& a, const double& b)
    : a(a), b(b) {/* empty */}

// compute the force
void DoubleWell::compute_force(Molecule* molecule_pt) {

    // zero the potential
    molecule_pt->potential() = 0.0;
    // zero the laplacian
    molecule_pt->laplace() = 0.0;
    // zero the magnetisation
    molecule_pt->magnetisation() = 0.0;

    for(auto& particle : molecule_pt->Particles) {
        // calculate the new force
        particle.second->f = get_force(particle.second->q);

        // add contribution to potential
        molecule_pt->potential() += get_potential(particle.second->q);
        // add contribution to laplacian
        molecule_pt->laplace() += get_laplace(particle.second->q);
        // add contribution to magnetisation
        molecule_pt->magnetisation() += get_applied_field(particle.second->q);
    }
}

// compute the force
Vector DoubleWell::compute_force(Molecule* molecule_pt,
                     Particle* particle_i,
                     Particle* particle_j,
                     const double& r,
                     const Vector& dr) {
                         // empty as not needed
                     }

// function for the potential
double DoubleWell::get_potential(const double& x) {
    return a * pow(pow(x, 2.0) - 1.0, 2.0) + b * (x + 1.0);
}

// function which returns the force Vector
double DoubleWell::get_potential(const Vector& x) {
    // make scalar to return
    double V = 0.0;

    // calculate the new force
    for(unsigned i=0; i<x.size(); i++) {
        V += get_potential(x(i));
    }

    // return the potential scalar
    return V;
}

// function which returns the applied field
double DoubleWell::get_applied_field(const double& x) {
    return x + 1.0;
}

// function which retruns the applied field
double DoubleWell::get_applied_field(const Vector& x) {
    // make scalar to return
    double B = 0.0;

    for(unsigned i=0; i<x.size(); i++) {
        B += get_applied_field(x(i));
    }

    return B;
}

// function for the laplacian
double DoubleWell::get_laplace(const double& x) {
    return 4.0 * a * (3.0 * pow(x, 2.0) - 1);
}

// function which returns the force Vector
double DoubleWell::get_laplace(const Vector& x) {
    // make scalar to return
    double L = 0.0;

    // calculate the new force
    for(unsigned i=0; i<x.size(); i++) {
        L += get_laplace(x(i));
    }

    // return the potential scalar
    return L;
}

// function which returns the force scalar
double DoubleWell::get_force(const double& x) {
    return -4.0 * a * x * (pow(x, 2.0) - 1) - b;
}

// function which returns the force Vector
Vector DoubleWell::get_force(const Vector& x) {
    // make vector to return
    Vector F = Vector(x.size());

    // calculate the new force
    for(unsigned i=0; i<x.size(); i++) {
        F(i) = get_force(x(i));
    }

    // return the force vector
    return F;
}
