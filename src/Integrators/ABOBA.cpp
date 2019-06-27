#include "ABOBA.hpp"

using namespace::std;

/******************************************************************************
                                 ABOBA Class
 *****************************************************************************/
ABOBA::ABOBA(const double& beta, const double& gamma, const double& gamma_rot,
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
ABOBA::~ABOBA()
{
    delete &Time_Step;
    delete &Step;

    delete &With_npt;
    delete &With_st;
}

// integrator
void ABOBA::integrate(Molecule* molecule_pt)
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
void ABOBA::integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
                                    Molecule* molecule_pt, const double& mass,
                                    const double& target_press,
                                    const double& gamma_box, const int& recf,
                                    const int& rect)
{
    With_npt = true;

    Langevin::integrate_with_npt_grid(Szero, cut_off, molecule_pt,
                                      mass, target_press, gamma_box, Time_Step,
                                      recf, rect);
}

// set with simulation tempering
void ABOBA::integrate_with_st(const double& tmin,
                              const double& tmax,
                              const double& n_temperatures,
                              const unsigned& mod_switch,
                              const int& seed)
{
    With_st = true;

    // initialise the simulated tempering class
    Langevin::integrate_with_st(tmin, tmax, n_temperatures, mod_switch, seed);
}

void ABOBA::set_npt_integrator_version(const unsigned& version)
{
    Npt_version = version;
}

// Integrate forward using NVT
void ABOBA::nvt_integration(Molecule* molecule_pt)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();
    Particle* particle = NULL;

    // preforce integration
    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        Langevin::A(*particle, 0.5 * Time_Step);
        Langevin::B(*particle, 0.5 * Time_Step);
        Langevin::O(*particle);
    }

    // force solve
    Langevin::compute_force(molecule_pt);

    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle = &molecule_pt->particle(i);

        // integrate forward
        Langevin::B(*particle, 0.5 * Time_Step);
        Langevin::A(*particle, 0.5 * Time_Step);
    }
}

// Integrate forward using NPT
void ABOBA::npt_integration(Molecule* molecule_pt)
{
    double before = 0.0;
    double after = 0.0;

    if(Npt_version == 1)
    {
        // integrate forward
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::O_NPT(molecule_pt);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);

    }
    else if(Npt_version == 2)
    {
        // integrate forward
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::O_NPT(molecule_pt);

        // force solve
        Langevin::compute_force(molecule_pt);

        // integrate forward
        Langevin::B_NPT(molecule_pt, 0.5 * Time_Step);
        Langevin::A_1_NPT(0.5 * Time_Step);
        Langevin::A_2_NPT(molecule_pt, 0.5 * Time_Step);
    }

}
