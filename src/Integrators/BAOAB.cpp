#include "BAOAB.hpp"

using namespace::std;

/******************************************************************************
                                 BAOAB Class
 *****************************************************************************/
BAOAB::BAOAB(const double& beta, const double& gamma, const double& gamma_rot,
             const double& time_step, System* system_pt, const int& seed)
             : Langevin(beta, gamma, gamma_rot, time_step, system_pt, seed)
{
    Time_Step = time_step;
    Step = 0;
    Npt_version = 1;

    With_npt = false;
    With_st = false;
}

// destructor
BAOAB::~BAOAB()
{
    delete &Time_Step;
    delete &Step;

    delete &With_npt;
    delete &With_st;
}

// integrator
void BAOAB::integrate(Molecule* molecule_pt)
{
    // update the step counter
    Step++;

    if(With_npt)
        npt_integration(molecule_pt);
    else
        nvt_integration(molecule_pt);

    if(With_st)
        Langevin::update_simulated_tempering(molecule_pt, Step, Time_Step);
}

// set with npt grid
void BAOAB::integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
                                    Molecule* molecule_pt, const double& mass,
                                    const double& target_press,
                                    const double& gamma_box)
{
    With_npt = true;

    Langevin::integrate_with_npt_grid(Szero, cut_off, molecule_pt,
                                      mass, target_press, gamma_box, Time_Step);
}

// set with simulation tempering
void BAOAB::integrate_with_st(const double& tmin,
                              const double& tmax,
                              const double& n_temperatures,
                              const unsigned& mod_switch,
                              const int& seed)
{
    With_st = true;

    // initialise the simulated tempering class
    Langevin::integrate_with_st(tmin, tmax, n_temperatures, mod_switch, seed);
}

void BAOAB::set_npt_integrator_version(const unsigned& version)
{
    Npt_version = version;
}

// Integrate forward using NVT
void BAOAB::nvt_integration(Molecule* molecule_pt)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();
    Particle* particle = NULL;

    // preforce integration
    // #pragma omp barrier
    // #pragma omp single
    // #pragma omp parallel for simd
    // #pragma omp parallel for default(shared) firstprivate(number_of_particles)
    for(unsigned i=0; i<number_of_particles; i++)
    // for(auto& particle = molecule_pt->Particles.begin();
    //     particle != molecule_pt->Particles.end(); particle++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        // printf("B step 1\n");
        Langevin::B(*particle, 0.5 * Time_Step);
        // printf("A step 1\n");
        Langevin::A(*particle, 0.5 * Time_Step);
        // printf("O step 1\n");
        Langevin::O(*particle);
        // printf("A step 2\n");
        Langevin::A(*particle, 0.5 * Time_Step);
    }

    // force solve
    // printf("\tForce\n");
    Langevin::compute_force(molecule_pt);

    // #pragma omp barrier
    // #pragma omp single
    // #pragma omp parallel for simd
    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        // printf("B step 2\n");
        Langevin::B(*particle, 0.5 * Time_Step);
    }
}

// Integrate forward using NPT
void BAOAB::npt_integration(Molecule* molecule_pt)
{
    double before = 0.0;
    double after = 0.0;

    if(Npt_version == 1)
    {
        // integrate forward
        // before = omp_get_wtime();
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        // after = omp_get_wtime();
        // printf("B step %1.5f\n", after-before);

        // before = omp_get_wtime();
        Langevin::A_1_NPT(0.5 * Time_Step);
        // after = omp_get_wtime();
        // printf("A1 step %1.5f\n", after-before);

        // before = omp_get_wtime();
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        // after = omp_get_wtime();
        // printf("A2 step %1.5f\n", after-before);

        // before = omp_get_wtime();
        Langevin::O_NPT(molecule_pt);
        // after = omp_get_wtime();
        // printf("O step %1.5f\n", after-before);

        // before = omp_get_wtime();
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        // after = omp_get_wtime();
        // printf("A step %1.5f\n", after-before);

        // force solve
        // before = omp_get_wtime();
        Langevin::compute_force(molecule_pt);
        // after = omp_get_wtime();
        // printf("Force step %1.5f\n", after-before);

        // integrate forward
        // before = omp_get_wtime();
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        // after = omp_get_wtime();
        // printf("B step %1.5f\n", after-before);

        // for(unsigned i=0; i<molecule_pt->nparticle(); i++)
        // {
        //     printf("\tParticle %d q %f %f\n",i,molecule_pt->particle(0).q(0),molecule_pt->particle(0).q(1) );
        //     printf("\tParticle %d p %f %f\n",i,molecule_pt->particle(0).p(0),molecule_pt->particle(0).p(1) );
        //     printf("\tParticle %d f %f %f\n",i,molecule_pt->particle(0).f(0),molecule_pt->particle(0).f(1) );
        // }
    }
    else if(Npt_version == 2)
    {
        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::O_NPT(molecule_pt);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
    }

}
