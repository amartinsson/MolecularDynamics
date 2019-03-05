#include "HarmonicOscillator.hpp"

/******************************************************************************
                           Harmonic Oscillator
 *****************************************************************************/
// Constructor
HarmonicOscillator::HarmonicOscillator(Molecule* molecule_pt)
{
    compute_force(molecule_pt);
};

// compute the force
void HarmonicOscillator::compute_force(Molecule* molecule_pt)
{
    // dereference helpers
    double DIM = molecule_pt->dim();
    unsigned nparticles = molecule_pt->nparticle();
    Particle* particle_i = NULL;

    double V = 0.0;

#pragma omp parallel for schedule(static) reduction(+:V)
    for(unsigned i=0; i<nparticles; i++)
    {
        particle_i = &molecule_pt->particle(i);
        particle_i->f.zero();

        particle_i->f -= particle_i->q;
        V += particle_i->q.dot(particle_i->q);
    }

    molecule_pt->potential() = V;
}

// compute the pair force
vector<double> HarmonicOscillator::compute_pair_force(Molecule* molecule_pt,
                                                      Particle* particle_i,
                                                      Particle* particle_j,
                                                      const double& r,
                                                      const double& r_x,
                                                      const double& r_y)
{
    // dereference helpers
    vector<double> forces(2, 0.0);

    forces[0] = -r_x;
    forces[1] = -r_y;

    // calculate the force
#pragma omp atomic
    particle_i->f(0) -= forces[0];
#pragma omp atomic
    particle_i->f(1) -= forces[1];

#pragma omp atomic
    particle_j->f(0) += forces[0];
#pragma omp atomic
    particle_j->f(1) += forces[1];

    // assign the potential
#pragma omp atomic
    molecule_pt->potential() += r_x*r_x + r_y*r_y;

    return forces;
}
