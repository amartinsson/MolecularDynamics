#include "Langevin.hpp"

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
    with_is = false;
    with_dis = false;

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

    delete &With_grid;
    delete &With_npt;

    delete System_pt;
    delete St_pt;
}

// Langevin based position step
void Langevin::A(Particle& particle, const double& h)
{
    particle.q += particle.m.inv() * particle.p * h;

    if(particle.rigid_body())
        A_rot(particle, h);
}

// Langevin based position step, constant pressure
void Langevin::A_1_NPT(const double& h)
{
    // make references
    Matrix Sold = NptGrid_pt->S;

    double npt_mass = NptGrid_pt->get_mass();

    NptGrid_pt->S += NptGrid_pt->Sp * (h / npt_mass);

    // need to rescale the particles such that they are
    // in the same position with repsect to the cell
    NptGrid_pt->enforce_constant_relative_particle_pos(Sold);
}

// Langevin based position step, constant pressure
void Langevin::A_2_NPT(Molecule* molecule_pt, const double& h)
{
    // particle integration
    #pragma omp parallel for default(shared) schedule(static)
    for(unsigned i=0; i<molecule_pt->nparticle(); i++)
        A(molecule_pt->particle(i), h);

    // update box momentum
    A_NPT_2_box(h);
}


// Langevin based position step, constant pressure, box
void Langevin::A_NPT_2_box(const double& h)
{
    // update kinetic gradient
    NptGrid_pt->update_kinetic_gradient();

    // update the box momentum
    NptGrid_pt->Sp -= NptGrid_pt->nablaK * h;
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
    //#pragma omp simd
    #pragma omp parallel for default(shared) schedule(static)
    for(unsigned i=0; i<number_of_particles; i++)
        B(molecule_pt->particle(i), h);

    // box integration
    B_NPT_box(h);
}

void Langevin::B_NPT_box(const double& h)
{
    // make diagonal matrix
    Matrix Sd = NptGrid_pt->S;
    double det = NptGrid_pt->S.det();

    double Target_Pressure = NptGrid_pt->get_target_pressure();

    Sd = Sd.inv().off_diag_zero();

    NptGrid_pt->Sp -= (NptGrid_pt->virial + Sd * det * Target_Pressure) * h;
}

NptGrid& Langevin::npt_obj() {
    return *NptGrid_pt;
}

Grid& Langevin::grid_obj() {
    return *Grid_pt;
}

InfiniteSwitch& Langevin::is_obj() {
    return *is_pt;
}

DoubleInfiniteSwitch& Langevin::dis_obj() {
    return *dis_pt;
}

SimulatedTempering& Langevin::st_obj() {
    return *St_pt;
}

OngulatedTempering& Langevin::ot_obj() {
    return *Ot_pt;
}

SimulatedAnnealing& Langevin::sa_obj() {
    return *Sa_pt;
}

ReplicaExchange& Langevin::re_obj() {
    return *Re_pt;
}

void Langevin::npt_set_initial(Molecule* molecule_pt,
                               const char* initial_pos_filename,
                               const char* initial_mom_filename,
                               const char* initial_box_filename)
{
    if(With_npt)
    {
        // read in the particles and the box
        NptGrid_pt->add_file_initial_condition(molecule_pt,
                                               initial_pos_filename,
                                               initial_mom_filename,
                                               initial_box_filename);

        // check the initial condition
        NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
    }
    else
    {
        // read in the particles and the box
        Grid_pt->add_file_initial_condition(molecule_pt,
                                               initial_pos_filename,
                                               initial_mom_filename,
                                               initial_box_filename);

        // check the initial condition
        Grid_pt->update_particle_forces(System_pt, molecule_pt);
    }

}

// compute the force in the correct way
void Langevin::compute_force(Molecule* molecule_pt)
{

    if(With_grid) {
        Grid_pt->update_particle_forces(System_pt, molecule_pt);
    }
    else if(With_npt) {
        NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
    }
    else {
        System_pt->compute_force(molecule_pt);
    }

    // if we are running with ISST rescale the forces
    if(with_is) {
        is_pt->apply_force_rescaling();
    }

    if(with_dis) {
        dis_pt->apply_force_rescaling();
    }

    // // update the force for all the replicas
    // if(With_re) {
    //     unsigned N = Re_pt->replica.size();
    //
    //     for(unsigned i=1; i<N; i++) {
    //         System_pt->compute_force(Re_pt->replica(i));
    //     }
    // }
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

    // particle integration
    // WARNING cannot parallelisation on this step unless random
    // numbers are generated individially for each open mp thread
    #pragma omp simd
    for(unsigned i=0; i<number_of_particles; i++)
        O(molecule_pt->particle(i));

    // box integration
    O_NPT_box();
}

void Langevin::O_NPT_box()
{
    // get the dimension of momentum
    int dim = NptGrid_pt->Sp.size()[0];
    double npt_mass = NptGrid_pt->get_mass();

    // make an upper triangular random matrix
    Matrix N(dim, dim, normal_gen);

    NptGrid_pt->Sp = NptGrid_pt->Sp * OstepCbox
                    + N * OstepZbox * sqrt(npt_mass);
}

// Integrator function
void Langevin::integrate(Molecule* molecule_pt)
{
    std::cout << "Error: Langevin class must implement integrate() function\n";
    exit(-1);
}

// set with grid
void Langevin::integrate_with_grid(const Matrix& Szero, const double& cut_off,
                                   Molecule* molecule_pt, const int& recf,
                                   const int& rect)
{
    // set boolean
    With_grid = true;
    // make a new grid
    Grid_pt = new Grid(Szero, cut_off, molecule_pt, recf, rect);
    // update all the particle forces using this grid
    Grid_pt->update_particle_forces(System_pt, molecule_pt);
}

// set with npt_grid
void Langevin::integrate_with_npt_grid(const Matrix& Szero,
                                       const double& cut_off,
                                       Molecule* molecule_pt,
                                       const double& mass,
                                       const double& target_press,
                                       const double& gamma_npt,
                                       const double& o_box_time_step,
                                       const int& recf, const int& rect)
{
    // set parameters
    With_npt = true;
    // Target_pressure = target_press;
    GammaNpt = gamma_npt;
    // make new npt grid
    NptGrid_pt = new NptGrid(Szero, cut_off, molecule_pt, mass,
                             target_press, recf, rect);
    // update all the particle forces on this grid
    NptGrid_pt->update_particle_forces(System_pt, molecule_pt);
    // update the box forces
    // NptGrid_pt->update_box_force();

    // Box o step values
    OstepCbox = exp(-GammaNpt * o_box_time_step);
    OstepZbox = sqrt((1.0 - OstepCbox * OstepCbox) * 1.0 / Beta);
}


void Langevin::integrate_with_infinite_switch(InfiniteSwitch* scheme)
{
    with_is = true;

    is_pt = scheme;

    // if we're integrating with npt, initialize isst with it.
    if(With_npt) {
        is_pt->initialize_with_npt(NptGrid_pt);
    }
}

void Langevin::integrate_with_double_infinite_switch(DoubleInfiniteSwitch* scheme)
{
    with_dis = true;

    dis_pt = scheme;
}

// set with simulation tempering
void Langevin::integrate_with_st(const double& tmin,
                                 const double& tmax,
                                 const double& n_temperatures,
                                 const unsigned& mod_switch,
                                 const int& seed)
{
    With_st = true;
    // initialise the simulated tempering class
    St_pt = new SimulatedTempering(tmin, tmax, n_temperatures,
                                   mod_switch, seed);
}

// set with ongulated tempering
void Langevin::integrate_with_ot(const double& tmin,
                                 const double& tmax,
                                 const double& n_temperatures,
                                 const unsigned& mod_switch,
                                 const int& seed)
{
    With_ot = true;
    // initialise the simulated tempering class
    Ot_pt = new OngulatedTempering(tmin, tmax, n_temperatures,
                                   mod_switch, seed);
}

// set with Simulated Annealing
void Langevin::integrate_with_sa(const double& tmax, const double& crate,
    const unsigned& nsteps, Molecule* molecule_pt)
{
    With_sa = true;
    // initialise the class
    Sa_pt = new SimulatedAnnealing(tmax, crate, nsteps, molecule_pt);
}

void Langevin::integrate_with_re(const double& tmin, const double& tmax,
    Molecule* molecule_pt, const unsigned& Ncopies, const unsigned& sfreq,
        const int& seed)
{
    With_re = true;
    // initialise the class
    Re_pt = new ReplicaExchange(tmin, tmax, molecule_pt, Ncopies, sfreq, seed);
}

// check if to update temperature
void Langevin::update_simulated_tempering(Molecule* molecule_pt,
                                          const unsigned& step,
                                          const double& h)
{
    // update hte temperature and all the variables
    St_pt->update_temperature(molecule_pt, step);
    update_integrator_temperature(1.0 / St_pt->get_temperature(), h);
}

// check if to update temperature
void Langevin::update_ongulated_tempering(Molecule* molecule_pt,
                                          const unsigned& step,
                                          const double& h)
{
    // update all the temperatures and variables
    Ot_pt->update_temperature(molecule_pt, step);
    update_integrator_temperature(1.0 / Ot_pt->get_temperature(), h);
}

void Langevin::update_simulated_annealing(const unsigned& step,
    const double& h)
{
    Sa_pt->update_temperature(step);
    update_integrator_temperature(Sa_pt->get_beta(), h);
}

void Langevin::update_replica_exchange(const double& h)
{
    printf("WARNGIN: Rep. Exch. not implimented for method!\n");
}

// update temperature
void Langevin::update_integrator_temperature(const double& beta,
        const double& h) {

    // set beta
    Beta = beta;

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

    // calculate the new angle
    pi = particle.pi(0,0);
    I = particle.I(0,0);
    alpha = h * pi / I;

    // make rotation matrix
    RotMatrix R(alpha);

    // multiply matrices together
    Matrix Q_np1 = particle.Q * R.T();

    // update the Q matrix
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

// check if cell blow up
bool Langevin::cell_blow_up()
{
    bool state = false;

    if(With_npt)
        if(state != NptGrid_pt->break_experiment)
            state = NptGrid_pt->break_experiment;

    // return the state
    return state;
}
