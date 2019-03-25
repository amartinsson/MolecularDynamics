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
void Langevin::A(Particle& particle, const double& h)
{
    particle.q += (particle.p.T() * particle.m.inv()) * h;

    if(particle.rigid_body())
        A_rot(particle, h);
}

// Langevin based position step, constant pressure
void Langevin::A_1_NPT(const double& h)
{
    // make references
    Matrix Sold = NptGrid_pt->S;

    NptGrid_pt->S += NptGrid_pt->Sp * (h / Npt_mass);

    // need to rescale the particles such that they are
    // in the same position with repsect to the cell
    NptGrid_pt->enforce_constant_relative_particle_pos(Sold);

    // update the gradient matrices
    NptGrid_pt->update_gradient_matrices();
}

// Langevin based position step, constant pressure
void Langevin::A_2_NPT(Molecule* molecule_pt, const double& h)
{
    // particle integration
    #pragma omp simd
    for(unsigned i=0; i<molecule_pt->nparticle(); i++)
    {
        // position dependent part
        A_NPT_2_part(molecule_pt->particle(i), h);

        // particle add contribution to box
        A_NPT_2_box(molecule_pt->particle(i), h);
    }
}

// Langevin based position step, constant pressure, partilce
void Langevin::A_NPT_2_part(Particle& particle, const double& h)
{
    // update particle position
    particle.q += NptGrid_pt->S * NptGrid_pt->S.inv().T() * particle.m.inv()
                * particle.p * h;

    if(particle.rigid_body())
        A_rot(particle, h);
}

// Langevin based position step, constant pressure, box
void Langevin::A_NPT_2_box(Particle& particle, const double& h)
{
    // copy matrices
    Matrix Ka = NptGrid_pt->nablaSa / particle.m(0,0);
    Matrix Kb = NptGrid_pt->nablaSb / particle.m(0,0);
    Matrix Kd = NptGrid_pt->nablaSd;

    // divide by the correct mass
    Kd(0,1) /= particle.m(0,0);
    Kd(1,0) /= particle.m(0,0);
    Kd(1,1) /= particle.m(1,1);

    // update the box momentum correctly
    NptGrid_pt->Sp(0,0) -= h * (particle.p.T() * Ka).dot(particle.p);
    NptGrid_pt->Sp(0,1) -= h * (particle.p.T() * Kb).dot(particle.p);
    NptGrid_pt->Sp(1,1) -= h * (particle.p.T() * Kd).dot(particle.p);
}

// Langevin Momentum based update
void Langevin::B(Particle& particle, const double& h)
{
    // update momentum
    particle.p += particle.f * h;

    if(particle.rigid_body())
        B_rot(particle, h);
}

// Langevin Momentum based update, constant pressure
void Langevin::B_NPT(Molecule* molecule_pt, const double& h)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();

    // particle integration
    #pragma omp simd
    for(unsigned i=0; i<number_of_particles; i++)
        B_NPT_part(molecule_pt->particle(i), h);

    // box integration
    B_NPT_box(h);
}

void Langevin::B_NPT_part(Particle& particle, const double& h)
{
    particle.p += NptGrid_pt->S.inv() * NptGrid_pt->S.T() * particle.f * h;

    if(particle.rigid_body())
        B_rot(particle, h);
}


void Langevin::B_NPT_box(const double& h)
{
    // make diagonal matrix
    Matrix Sd = NptGrid_pt->S;
    double det = NptGrid_pt->S.det();

    Sd = Sd.inv();
    Sd(0,1) = 0.0;

    NptGrid_pt->Sp -= (NptGrid_pt->virial + Sd * det * Target_pressure) * h;
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
void Langevin::O(Particle& particle)
{
    // generate random numbers
    Vector N(particle.p.size(), normal_gen);

    // find new momentum
    particle.p = particle.p * OstepC + particle.m.sqrt() * N * OstepZ;

    if(particle.rigid_body())
        O_rot(particle);
}

// Langevin Stochastic Momentum based update, constant pressure
void Langevin::O_NPT(Molecule* molecule_pt)
{
    // helper dereference
    unsigned number_of_particles = molecule_pt->nparticle();

    // WARNING
    // Warning this is bad and assumes that all the particles
    // have the same mass!!!!
    // WARNING

    // varying mass
    Matrix MassMatrix = NptGrid_pt->S.inv()
                * (NptGrid_pt->S * molecule_pt->particle(0).m
                    * NptGrid_pt->S.T()).symsqrt();

    // particle integration
    #pragma omp simd
    for(unsigned i=0; i<number_of_particles; i++)
        O_NPT_part(molecule_pt->particle(i), MassMatrix);

    // box integration
    O_NPT_box();
}

void Langevin::O_NPT_part(Particle& particle, Matrix& MassMatrix)
{
    // generate random numbers
    Vector N(particle.p.size(), normal_gen);

    // find new momentum
    particle.p = particle.p * OstepC + MassMatrix * N * OstepZ;

    if(particle.rigid_body())
        O_rot(particle);
}

void Langevin::O_NPT_box()
{
    // make an upper triangular random matrix
    Matrix N(2, 2, normal_gen);

    NptGrid_pt->Sp = NptGrid_pt->Sp * OstepCbox
                    + N * OstepZbox * sqrt(Npt_mass);
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

void Langevin::A_rot(Particle &particle, const double& h)
{
    // dereference helpers
    unsigned DIM = particle.dim();

    if(DIM > 2)
    {
        std::cout << "Error: Rotation can only be used in 2D\n";
        exit(-1);
    }

    // dereference helpers
    double pi = 0.0;
    double I = 0.0;
    double alpha = 0.0;

    pi = particle.pi(0,0);
    I = particle.I(0,0);
    alpha = h * pi / I;

    printf("ERROR: NEED TO CHECK THAT THIS IS CORRECT -- ARE THEY SPINNING IN THE RIGHT DIRECTION\n");
    exit(-1);

    RotMatrix R(alpha);

    // clockwise spinning:
    // R^T * Q
    // anti-clockwise spinning:
    // R * Q

    // multiply matrices together
    Matrix Q_np1 = R.T() * particle.Q;
    particle.Q = Q_np1;
}

void Langevin::B_rot(Particle& particle, const double& h)
{
    // update momentum
    particle.pi += particle.tau * h;
}

void Langevin::O_rot(Particle& particle)
{
    particle.pi(0,0) = OstepCphi * particle.pi(0,0)
                     + OstepZphi * sqrt(particle.I(0,0)) * normal_gen();
}
