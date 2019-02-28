#include "Integrators.hpp"

using namespace::std;

/******************************************************************************
                           Langevin Base Class
 *****************************************************************************/
Langevin::Langevin(const double& beta, const double& gamma,
                   const double& gamma_rot, const double& o_time_step,
                   System* system_pt, const int& seed)
                   : normal_gen(0.0, 1.0, seed)
{
    // parameters
    Beta = beta;
    Gamma = gamma;
    GammaRot = gamma_rot;

    // set force calculation
    System_pt = system_pt;

    // set booleans
    With_isst = false;
    With_grid = false;
    With_npt = false;

    // translational
    OstepC = exp(-Gamma * o_time_step);
    OstepZ = sqrt((1.0 - OstepC * OstepC) * 1.0 / Beta);

    // rotational
    OstepCphi = exp(-GammaRot * o_time_step);
    OstepZphi = sqrt((1.0 - OstepCphi * OstepCphi) * 1.0 / Beta);
}

// destructor
Langevin::~Langevin()
{
    delete[] &Beta;
    delete[] &Gamma;
    delete[] &GammaRot;

    delete[] &OstepC;
    delete[] &OstepZ;
    delete[] &OstepCphi;
    delete[] &OstepZphi;
    delete[] &OstepCbox;
    delete[] &OstepZbox;

    delete &Target_pressure;
    delete &Npt_mass;

    delete &With_isst;
    delete &With_grid;
    delete &With_npt;

    delete System_pt;
    delete St_pt;
}

// Langevin based position step
void Langevin::A(Particle* particle_pt, const double& h)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();
    double* q = NULL;
    double  p = 0.0;
    double  m = 0.0;

    for(unsigned j=0; j<DIM; j++)
    {
        q = particle_pt->q_pt(j);
	    p = *particle_pt->p_pt(j);
	    m = *particle_pt->m_pt(j);

        *q += h * p / m;
    }

    if(particle_pt->rigid_body() == true)
        A_rot(particle_pt, h);
}

// Langevin based position step, constant pressure
void Langevin::A_1_NPT(const double& h)
{
    // make references
    double* L = NULL;
    double Lp = 0.0;
    double mass = NptGrid_pt->get_mass();
    unsigned number_of_grid_coordinates = NptGrid_pt->get_ncoord();

    // old position needed for rescaling
    vector<double> L_old(3, 0.0);

    for(unsigned i=0; i<number_of_grid_coordinates; i++)
    {
        // derefernce the grid coordinates
        L = NptGrid_pt->L_pt(i);
        Lp = *NptGrid_pt->Lp_pt(i);
        // save the old position
        L_old[i] = *L;
        // update the coorinate
        *L += h * Lp / mass;
    }

    // need to rescale the particles such that they are
    // in the same position with repsect to the cell
    NptGrid_pt->enforce_constant_relative_particle_pos(L_old);
}

// Langevin based position step, constant pressure
void Langevin::A_2_NPT(Molecule* molecule_pt, const double& h)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();
    Particle* particle_pt = NULL;

    // particle integration
    #pragma omp parallel for simd
    for(unsigned i=0; i<number_of_particles; i++)
    {
        particle_pt = molecule_pt->particle_pt(i);
        A_NPT_2_part(particle_pt, h);
    }

    // box integration
    A_NPT_2_box(h);
}

// Langevin based position step, constant pressure, partilce
void Langevin::A_NPT_2_part(Particle* particle_pt, const double& h)
{
    // holder for particle positon
    double* q = NULL;
    unsigned DIM = particle_pt->dim();

    // dereference the box size
    double lx1 = *NptGrid_pt->L_pt(0);
    double ly1 = *NptGrid_pt->L_pt(1);
    double ly2 = *NptGrid_pt->L_pt(2);

    // get the momentum of the particle
    double px = *particle_pt->p_pt(0);
    double py = *particle_pt->p_pt(1);

    // get the momentum of the particle
    double mx = *particle_pt->m_pt(0);
    double my = *particle_pt->m_pt(1);

    // loop over the number of dimensions
    for(unsigned j=0; j<DIM; j++)
    {
        // dereference the particle position
        q = particle_pt->q_pt(j);

        if(j==0)
        {
            (*q) += h *((1.0 - ly1*ly1/(lx1*ly2)) * px/mx
            + ly1/ly2 * py/my);
        }
        else if(j==1)
        {
            (*q) += h *(py/my - ly1/lx1 * px/mx);
        }
    }

    if(particle_pt->rigid_body() == true)
        A_rot(particle_pt, h);
}

// Langevin based position step, constant pressure, box
void Langevin::A_NPT_2_box(const double& h)
{
    // dereference the box size
    double lx1 = *NptGrid_pt->L_pt(0);
    double ly1 = *NptGrid_pt->L_pt(1);
    double ly2 = *NptGrid_pt->L_pt(2);

    // update the accumulated momentum in the box
    NptGrid_pt->update_accumulted_momentum();
    double acc_mom = 0.0;
    double* Lp = NULL;
    unsigned number_of_grid_coordinates = NptGrid_pt->get_ncoord();

    // loop over all the box coordinates
    for(unsigned i=0; i<number_of_grid_coordinates;i++)
    {
        acc_mom = NptGrid_pt->get_accumulated_momentum(i);
        Lp = NptGrid_pt->Lp_pt(i);

        if(i != 2)
            (*Lp) += h * (acc_mom / lx1);
        else
        {
            double acc_mom_xy = NptGrid_pt->get_accumulated_momentum(1);
            (*Lp) += h * (acc_mom / ly2 - ly1/(ly2*lx1) * acc_mom_xy);
        }
    }
}

// Langevin Momentum based update
void Langevin::B(Particle* particle_pt, const double& h)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();
    double* p = NULL;
    double f = 0.0;

    // loop over the number of dimensions
    for(unsigned j=0; j<DIM; j++)
    {
        p = particle_pt->p_pt(j);
	    f = *particle_pt->f_pt(j);

        *p += h * f;
    }

    if(particle_pt->rigid_body() == true)
        B_rot(particle_pt, h);
}

// Langevin Momentum based update, constant pressure
void Langevin::B_NPT(Molecule* molecule_pt, const double& h)
{
    // #pragma omp barrier
    // #pragma omp single
    {
        // helper dereference
        unsigned number_of_particles = molecule_pt->nparticle();
        Particle* particle_pt = NULL;

        // particle integration
        #pragma omp parallel for simd
        for(unsigned i=0; i<number_of_particles; i++)
        {
            particle_pt = molecule_pt->particle_pt(i);

            B_NPT_part(particle_pt, h);
        }

        // box integration
        B_NPT_box(h);
    }
}

void Langevin::B_NPT_part(Particle* particle_pt, const double& h)
{
    // holder for particle momentum
    double* p = NULL;
    unsigned DIM = particle_pt->dim();

    // dereference the box size
    double lx1 = *NptGrid_pt->L_pt(0);
    double ly1 = *NptGrid_pt->L_pt(1);
    double ly2 = *NptGrid_pt->L_pt(2);

    // get the forces on the particle
    double fx = *particle_pt->f_pt(0);
    double fy = *particle_pt->f_pt(1);

    for(unsigned j=0; j<DIM; j++)
    {
        // derefernce the particle momentum
        p = particle_pt->p_pt(j);

        if(j==0)
        {
            (*p) += h * ((1.0 - ly1*ly1/(lx1*ly2)) * fx - ly1/lx1 * fy);
        }
        else if(j==1)
        {
            (*p) += h * (fy + ly1/ly2 * fx);
        }
    }

    if(particle_pt->rigid_body() == true)
        B_rot(particle_pt, h);
}


void Langevin::B_NPT_box(const double& h)
{
    // dereference the box size
    double lx1 = *NptGrid_pt->L_pt(0);
    double ly1 = *NptGrid_pt->L_pt(1);
    double ly2 = *NptGrid_pt->L_pt(2);

    double* Lp = NULL;
    double virial = 0.0;
    unsigned number_of_grid_coordinates = NptGrid_pt->get_ncoord();

    for(unsigned i=0; i<number_of_grid_coordinates; i++)
    {
        Lp = NptGrid_pt->Lp_pt(i);
        virial = NptGrid_pt->get_virial(i);

        if(i == 0)
            *Lp -= h * (virial + Target_pressure * ly2);
        else if(i == 1)
            *Lp -= h * virial;
        else if(i == 2)
            *Lp -= h * (virial + Target_pressure * lx1);
    }
}

// update the pressure and temperature reading
void Langevin::npt_update_pressure_temperature()
{
    NptGrid_pt->update_pressure_temperature();
}

// get the instant pressure reading
double Langevin::npt_get_instant_pressure()
{
    return NptGrid_pt->get_instant_pressure();
}

// get the pressure reading
double Langevin::npt_get_pressure()
{
    return NptGrid_pt->get_pressure();
}

// get the instant temperature reading
double Langevin::npt_get_instant_temperature()
{
    return NptGrid_pt->get_instant_temperature();
}

// get the temperature reading
double Langevin::npt_get_temperature()
{
    return NptGrid_pt->get_temperature();
}

void Langevin::npt_set_initial(Molecule* molecule_pt,
                               const char* initial_pos_filename,
                               const char* initial_mom_filename,
                               const char* initial_box_filename)
{
    // read in the particles and the box
    NptGrid_pt->add_file_initial_condition(molecule_pt, initial_pos_filename,
        initial_mom_filename,
        initial_box_filename);

    // check the initial condition
    NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
}

// compute the force in the correct way
void Langevin::compute_force(Molecule* molecule_pt)
{

    if(With_grid)
    {
        Grid_pt->update_particle_forces(System_pt, molecule_pt);
    }
    else if(With_npt)
    {
        NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
    }
    else
    {
        System_pt->compute_force(molecule_pt);
    }

    // if we are running with ISST rescale the forces =
    if(With_isst)
        Isst_pt->apply_force_rescaling(molecule_pt);
}

// Langevin Stochastic Momentum based update
void Langevin::O(Particle* particle_pt)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();
    double* p = NULL;
    double m  = 0.0;

    // loop over the number of dimensions
    for(unsigned j=0; j<DIM; j++)
    {
        p = particle_pt->p_pt(j);
	    m = *particle_pt->m_pt(j);

        *p = OstepC * (*p) + OstepZ * sqrt(m) * normal_gen();
    }

    if(particle_pt->rigid_body() == true)
        O_rot(particle_pt);
}

// Langevin Stochastic Momentum based update, constant pressure
void Langevin::O_NPT(Molecule* molecule_pt)
{
    // #pragma omp barrier
    // #pragma omp single
    {
        // helper dereference
        unsigned number_of_particles = molecule_pt->nparticle();
        Particle* particle_pt = NULL;

        // particle integration
        #pragma omp parallel for simd
        for(unsigned i=0; i<number_of_particles; i++)
        {
            particle_pt = molecule_pt->particle_pt(i);

            O(particle_pt);
        }

        // box integration
        O_NPT_box();
    }
}

void Langevin::O_NPT_box()
{
    // dereference the box coordinates
    double* Lp;
    double box_mass = NptGrid_pt->get_mass();
    unsigned number_of_grid_coordinates = NptGrid_pt->get_ncoord();

    // loop over the box coordinates
    for(unsigned i=0; i<number_of_grid_coordinates; i++)
    {
        // dereference
        Lp = NptGrid_pt->Lp_pt(i);

        // update
        *Lp = OstepCbox * (*Lp) + OstepZbox * sqrt(box_mass) * normal_gen();
    }
}

// Integrator function
void Langevin::integrate(Molecule* molecule_pt)
{
    std::cout << "Error: Langevin class must implement integrate() function\n";
    exit(-1);
}

// set with grid
void Langevin::integrate_with_grid(const double& a_x, const double& b_x,
                                   const double& b_y, const double& cut_off,
                                   Molecule* molecule_pt)
{
    // set boolean
    With_grid = true;
    // make a new grid
    Grid_pt = new Grid(a_x, b_x, b_y, cut_off, molecule_pt);
    // update all the particle forces using this grid
    Grid_pt->update_particle_forces(System_pt, molecule_pt);
}

// set with npt_grid
void Langevin::integrate_with_npt_grid(const double& a_x, const double& b_x,
                                       const double& b_y, const double& cut_off,
                                       Molecule* molecule_pt,
                                       const double& mass,
                                       const double& target_press,
                                       const double& gamma_npt,
                                       const double& o_box_time_step)
{
    // set parameters
    With_npt = true;
    Target_pressure = target_press;
    Npt_mass = mass;
    GammaNpt = gamma_npt;
    // make new npt grid
    NptGrid_pt = new NptGrid(a_x, b_x, b_y, cut_off, molecule_pt, mass,
                             target_press);
    // update all the particle forces on this grid
    NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
    // update the box forces
    // NptGrid_pt->update_box_force();

    // Box o step values
    OstepCbox = exp(-GammaNpt * o_box_time_step);
    OstepZbox = sqrt((1.0 - OstepCbox * OstepCbox) * 1.0 / Beta);
}

// set with isst
void Langevin::integrate_with_isst(Molecule* molecule_pt, const double& tmin,
                                   const double& tmax, const unsigned& nint,
                                   const double& t_step, const double& tau)
{
    With_isst = true;
    Isst_pt = new InfiniteSwitchSimulatedTempering(molecule_pt, tmin, tmax,
                                                   nint, t_step, tau);
}

// set with simulation tempering
void Langevin::integrate_with_st(const double& tmin,
                                 const double& tmax,
                                 const double& n_temperatures,
                                 const unsigned& mod_switch,
                                 const int& seed)
{
    // initialise the simulated tempering class
    St_pt = new SimulatedTempering(tmin, tmax, n_temperatures,
                                   mod_switch, seed);
}

// check if to update temperature
void Langevin::update_simulated_tempering(Molecule* molecule_pt,
                                          const unsigned& step,
                                          const double& h)
{
    St_pt->update_temperature(molecule_pt, step);
    Beta = 1.0 / St_pt->get_temperature();

    // translational
    OstepC = exp(-Gamma * h);
    OstepZ = sqrt((1.0 - OstepC * OstepC) * 1.0 / Beta);

    // rotational
    OstepCphi = exp(-GammaRot * h);
    OstepZphi = sqrt((1.0 - OstepCphi * OstepCphi) * 1.0 / Beta);

    // npt
    if(With_npt)
    {
        OstepCbox = exp(-GammaNpt * h);
        OstepZbox = sqrt((1.0 - OstepCbox * OstepCbox) * 1.0 / Beta);
    }
}


void Langevin::A_rot(Particle* particle_pt, const double& h)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();

    if(DIM > 2)
    {
        std::cout << "Error: Rotation can only be used in 2D\n";
        exit(-1);
    }

    // dereference helpers
    double pi = 0.0;
    double I = 0.0;
    double alpha = 0.0;

    // rotation matrix
    double R_00 = 0.0;
    double R_01 = 0.0;
    double R_10 = 0.0;
    double R_11 = 0.0;

    for(unsigned i=0; i<1; i++)
    {
        pi = *particle_pt->pi_pt(i);
        I = *particle_pt->I_pt(i);
        alpha = h * pi / I;

        printf("ERROR: NEED TO CHECK THAT THIS IS CORRECT -- ARE THEY SPINNING IN THE RIGHT DIRECTION\n");
        exit(-1);

        RotMatrix R(alpha);

        // clockwise spinning:
        // R^T * Q
        // anti-clockwise spinning:
        // R * Q


        // multiply matrices together
        Matrix Q_np1 = R.T() * particle_pt->Q();
        particle_pt->Q() = Q_np1;

        // R = [  Cos(alpha) Sin(alpha)
        //       -Sin(alpha) Cos(alpha) ]
        // R_00 = cos(alpha);
        // R_01 = -sin(alpha);
        // R_10 = sin(alpha);
        // R_11 = cos(alpha);

        // Update the rotatin matrix
        //     Q = QR^T
        // double Q_00 = *particle_pt->Q_pt(0,0) * R_00
        //             + *particle_pt->Q_pt(0,1) * R_01;
        // double Q_01 = *particle_pt->Q_pt(0,0) * R_10
        //             + *particle_pt->Q_pt(0,1) * R_11;
        // double Q_10 = *particle_pt->Q_pt(1,0) * R_00
        //             + *particle_pt->Q_pt(1,1) * R_01;
        // double Q_11 = *particle_pt->Q_pt(1,0) * R_10
        //             + *particle_pt->Q_pt(1,1) * R_11;
        //
        // // Asign the new matrix
        // *particle_pt->Q_pt(0,0) = Q_00;
        // *particle_pt->Q_pt(0,1) = Q_01;
        // *particle_pt->Q_pt(1,0) = Q_10;
        // *particle_pt->Q_pt(1,1) = Q_11;
    }
}

void Langevin::B_rot(Particle* particle_pt, const double& h)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();

    if(DIM > 2)
    {
        std::cout << "Error: Rotation can only be used in 2D\n";
        exit(-1);
    }

    // dereference helpers
    double* pi = NULL;
    double tau = 0.0;

    for(unsigned j=0; j<1; j++)
    {
        pi = particle_pt->pi_pt(j);
        tau = *particle_pt->tau_pt(j);

        *pi += h * tau;
    }
}

void Langevin::O_rot(Particle* particle_pt)
{
    // dereference helpers
    unsigned DIM = particle_pt->dim();

    if(DIM > 2)
    {
        std::cout << "Error: Rotation can only be used in 2D\n";
        exit(-1);
    }

    // dereference helpers
    double* pi = NULL;
    double I = 0.0;

    for(unsigned j=0; j<DIM; j++)
    {
        pi = particle_pt->pi_pt(j);
        I = *particle_pt->I_pt(j);

        *pi = OstepCphi * (*pi) + OstepZphi * sqrt(I) * normal_gen();
    }
}
