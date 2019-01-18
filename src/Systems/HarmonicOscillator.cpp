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
    double* f = NULL;
    unsigned nparticles = molecule_pt->nparticle();
    double q = 0;
    Particle* particle_i = NULL;

    double* V = molecule_pt->potential_pt();
    *V = 0.0;

    #pragma omp for schedule(static)
    for(unsigned i=0; i<nparticles; i++)
    {
        particle_i = molecule_pt->particle_pt(i);

        for(unsigned j=0; j<DIM; j++)
        {
            f = particle_i->f_pt(j);
            q = *particle_i->q_pt(j);

            *f = -q;
            *V += q * q;
        }
    }
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
    double* V = molecule_pt->potential_pt();
    vector<double> forces(2, 0.0);

    forces[0] = -r_x;
    forces[1] = -r_y;

    // calculate the force
    #pragma omp atomic
    *particle_i->f_pt(0) -= forces[0];
    #pragma omp atomic
    *particle_i->f_pt(1) -= forces[1];

    #pragma omp atomic
    *particle_j->f_pt(0) += forces[0];
    #pragma omp atomic
    *particle_j->f_pt(1) += forces[1];

    // assign the potential
    #pragma omp atomic
    *V += r_x*r_x + r_y*r_y;

    return forces;
}
