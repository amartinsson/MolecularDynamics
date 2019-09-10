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
}

// destructor
BAOAB::~BAOAB()
{
    delete &Time_Step;
    delete &Step;
}

// set with npt grid
void BAOAB::integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
    Molecule* molecule_pt, const double& mass, const double& target_press,
        const double& gamma_box, const int& recf, const int& rect) {
            // call the Langevin grid integrator
            Langevin::integrate_with_npt_grid(Szero, cut_off,molecule_pt,
                mass, target_press, gamma_box, Time_Step, recf, rect);
}

// integrator
void BAOAB::integrate(Molecule* molecule_pt)
{
    // make sure integrator has correct temperature for replica exchange
    if(Langevin::With_re) {
        Langevin::update_integrator_temperature(
                            Langevin::Re_pt->replica(0)->beta(), Time_Step);
    }

    // update the step counter
    Step++;

    if(Langevin::With_npt)
        npt_integration(molecule_pt);
    else
        nvt_integration(molecule_pt);

    if(Langevin::With_st)
        Langevin::update_simulated_tempering(molecule_pt, Step, Time_Step);

    if(Langevin::With_ot)
        Langevin::update_ongulated_tempering(molecule_pt, Step, Time_Step);

    if(Langevin::With_sa)
        Langevin::update_simulated_annealing(Step, Time_Step);

    if(Langevin::With_re)
        this->update_replica_exchange(Time_Step);
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
    // int i = 32;
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
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);

        Langevin::A_1_NPT(0.25 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.25 * Time_Step);

        Langevin::O_NPT(molecule_pt);

        Langevin::A_1_NPT(0.25 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.25 * Time_Step);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
    }
    else if(Npt_version == 2)
    {
        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);

        Langevin::A_2_NPT(molecule_pt, 0.25 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.25 * Time_Step);

        Langevin::O_NPT(molecule_pt);

        Langevin::A_2_NPT(molecule_pt, 0.25 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.25 * Time_Step);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
    }

}

void BAOAB::update_replica_exchange(const double& h)
{
    unsigned N = Langevin::Re_pt->replica.size();

    // integrate forward for all replicas
    for(unsigned i=1; i<N; i++) {
        // get temperature
        double rep_beta = Langevin::Re_pt->replica(i)->beta();

        // set integrator temperature
        Langevin::update_integrator_temperature(rep_beta, h);

        // integrate replica forward
        this->nvt_integration(Langevin::Re_pt->replica(i));
    }

    // propose the moves
    Langevin::Re_pt->propose_moves(Step);

}
