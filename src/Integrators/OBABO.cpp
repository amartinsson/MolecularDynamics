#include "OBABO.hpp"

using namespace::std;

/******************************************************************************
                                 OBABO Class
 *****************************************************************************/
OBABO::OBABO(const double& beta, const double& gamma, const double& gamma_rot,
             const double& time_step, System* system_pt, const int& seed)
             : Langevin(beta, gamma, gamma_rot, 0.5 * time_step, system_pt, seed)
{
    Time_Step = time_step;
    Step = 0;
    Npt_version = 1;
}

// destructor
OBABO::~OBABO()
{
    delete &Time_Step;
    delete &Step;
}

// set with npt grid
void OBABO::integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
    Molecule* molecule_pt, const double& mass, const double& target_press,
        const double& gamma_box, const int& recf, const int& rect) {
            // call the Langevin grid integrator
            Langevin::integrate_with_npt_grid(Szero, cut_off,molecule_pt,
                mass, target_press, gamma_box, 0.5 * Time_Step, recf, rect);
}

// integrator
void OBABO::integrate(Molecule* molecule_pt)
{
    // update the step counter
    Step++;

    if(Langevin::With_npt)
        npt_integration(molecule_pt);
    else
        nvt_integration(molecule_pt);

    if(Langevin::With_st)
        Langevin::update_simulated_tempering(molecule_pt, 0.5 * Step, Time_Step);
}

void OBABO::set_npt_integrator_version(const unsigned& version)
{
    Npt_version = version;
}

// Integrate forward using NVT
void OBABO::nvt_integration(Molecule* molecule_pt)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();
    Particle* particle = NULL;

    // preforce integration
    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        Langevin::O(*particle);
        Langevin::B(*particle, 0.5 * Time_Step);
        Langevin::A(*particle, Time_Step);
    }

    // force solve
    Langevin::compute_force(molecule_pt);

    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        Langevin::B(*particle, 0.5 * Time_Step);
        Langevin::O(*particle);
    }
}

// Integrate forward using NPT
void OBABO::npt_integration(Molecule* molecule_pt)
{
    double before = 0.0;
    double after = 0.0;

    if(Npt_version == 1)
    {
        // integrate forward
        Langevin::O_NPT(molecule_pt);
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);

        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::O_NPT(molecule_pt);
    }
    else if(Npt_version == 2)
    {
        // integrate forward
        Langevin::O_NPT(molecule_pt);
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);

        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::O_NPT(molecule_pt);
    }
}
