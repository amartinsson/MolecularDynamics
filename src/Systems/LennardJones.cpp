#include "LennardJones.hpp"

using namespace::std;

LennardJones::LennardJones(const double& epsilon, const double& sigma) :
    Epsilon(epsilon), Sigma(sigma)
{
    // empty
}

// destructor
LennardJones::~LennardJones()
{
    delete[] &Epsilon;
    delete[] &Sigma;
}

// compute the force
void LennardJones::compute_force(Molecule* molecule_pt)
{
    std::cout << "Error: Lennard Jones force is a pair force\n";
    exit(-1);
}

Vector LennardJones::compute_force(Molecule* molecule_pt,
                                   Particle* particle_i,
                                   Particle* particle_j,
                                   const double& r,
                                   const Vector& dr)
{
    double dphidr = get_force_scalar(r);

    // we have removed minus coming from derivatie of x and y as the positive
    // direction is in going from i -> j and we need the force, which is the
    // negative of the gradient!
    Vector F = dr * dphidr / r;

    for(int k=0; k<F.size(); k++)
    {
        #pragma omp atomic
        particle_i->f(k) += F(k);

        #pragma omp atomic
        particle_j->f(k) -= F(k);
    }
    //
    //
    //
    // #pragma omp atomic
    // particle_i->f(1) += F(1);
    //
    // if(F.size() > 2)
    //     #pragma omp atomic
    //     particle_i->f(2) += F(2);
    //
    //
    //
    // #pragma omp atomic
    // particle_j->f(1) -= F(1);
    //
    // if(F.size() > 2)
    //     #pragma omp atomic
    //     particle_j->f(2) -= F(2);

    // add to the potential
    #pragma omp atomic
    molecule_pt->potential() += LennardJones::get_potential(r);

    return F;
}

// function for the Lennard Jones potential
double LennardJones::get_potential(const double& r)
{
    return 4.0 * Epsilon * (pow(Sigma / r, 12.0) - pow(Sigma / r, 6.0));
}

// function whcih returns the force scalar
double LennardJones::get_force_scalar(const double& r)
{
    return 24.0 * Epsilon * (pow(Sigma / r, 6.0) - 2.0 * pow(Sigma / r, 12.0)) / r;
}
