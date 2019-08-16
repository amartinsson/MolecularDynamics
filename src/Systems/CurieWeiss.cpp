#include "CurieWeiss.hpp"

using namespace::std;

CurieWeiss::CurieWeiss(const double& b, const unsigned& K)
    : b(b), K(K)
{
    // empyt
}

void CurieWeiss::compute_force(Molecule* molecule_pt)
{
    // calculate the magnetisation
    double m = compute_magnetisation(molecule_pt);
    // zero the potential
    molecule_pt->potential() = compute_potential(m);

    for(auto& particle : molecule_pt->Particles) {
        particle.second->f = get_force(m, particle.second->q);
    }
}

double CurieWeiss::compute_magnetisation(Molecule* molecule_pt)
{
    double m = 0.0;

    for(auto& particle : molecule_pt->Particles) {

        map_periodic(*particle.second);

        double theta = particle.second->q(0);

        m += cos(theta);
    }

    double magnetisation =  1.0/(double)K * m;

    molecule_pt->magnetisation() = magnetisation;

    return magnetisation;
}

double CurieWeiss::compute_potential(const double& m)
{
    return -(double)K * (0.5 * pow(m, 2.0) + b * m);
}

Vector CurieWeiss::get_force(const double& m, const Vector& theta)
{
    // only works for 1D so far
    Vector retf(1);

    retf(0) = -m * sin(theta(0)) - b * sin(theta(0));

    return retf;
}

void CurieWeiss::map_periodic(Particle& particle)
{
    if(fabs(particle.q(0)) >= 2.0 * M_PI) {
        particle.q(0) /= 2,0 * M_PI;
    }
}
