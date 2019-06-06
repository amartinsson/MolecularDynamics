#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

#include "Generator.hpp"
#include "Grid.hpp"
#include "InfiniteSwitchSimulatedTempering.hpp"
#include "NptGrid.hpp"
#include "Particle.hpp"
#include "SimulatedTempering.hpp"
#include "System.hpp"
#include "Array.hpp"

using namespace::std;

 /******************************************************************************
                            Generic Base Class
  *****************************************************************************/
class Integrator
{
public:
    Integrator() {};
    // destructor
    ~Integrator() {};
    // must be implemented
    virtual void integrate(Molecule* molecule_pt) = 0;
};

/******************************************************************************
                           Langevin Base Class
 *****************************************************************************/
class Langevin : public Integrator
{
public:
    Langevin(const double& beta, const double& gamma, const double& gamma_rot,
             const double& o_step_size, System* system_pt, const int& seed);
    // destructor
    ~Langevin();
    // integrator
    void integrate(Molecule* molecule_pt);
    // set to integrate with isst
    void integrate_with_isst(Molecule* molecule_pt, const double& tmin,
                             const double& tmax, const unsigned& nint,
                             const double& t_step, const double& tau);
    // set with grid
    void integrate_with_grid(const Matrix& Szero, const double& cut_off,
                             Molecule* molecule_pt, const int& recf,
                             const int& rect);
    // set the  initial condition from file
    void npt_set_initial(Molecule* molecule_pt,
                         const char* initial_pos_filename,
                         const char* initial_box_filename);
    // check blow up
    bool cell_blow_up();
    // return the npt grid observable
    NptGrid& npt_obj();
    Grid& grid_obj();

protected:
    void A(Particle& particle, const double& h);
    void B(Particle& particle, const double& h);
    void O(Particle& particle);

    void A_1_NPT(const double& h);
    void A_2_NPT(Molecule* molecule_pt, const double& h);
    void B_NPT(Molecule* molecule_pt, const double& h);
    void O_NPT(Molecule* molecule_pt);

    // compute the force in the correct way
    void compute_force(Molecule* molecule_pt);
    // set with npt grid
    void integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
                                 Molecule* molecule_pt, const double& mass,
                                 const double& target_press,
                                 const double& gamma_npt,
                                 const double& o_box_time_step,
                                 const int& recf, const int& rect);
    // set to integrate with simulated tempering
    void integrate_with_st(const double& tmin,
                           const double& tmax,
                           const double& n_temperatures,
                           const unsigned& mod_switch,
                           const int& seed);
    // check if to update temperature
    void update_simulated_tempering(Molecule* molecule_pt,
                                    const unsigned& step,
                                    const double& h);

private:
    double Beta;
    double Gamma;
    double GammaRot;
    double GammaNpt;

    double OstepC;
    double OstepZ;
    double OstepCphi;
    double OstepZphi;
    double OstepCbox;
    double OstepZbox;

    double Target_pressure;
    double Npt_mass;

    bool With_isst;
    bool With_grid;
    bool With_npt;

    NormalGenerator normal_gen;

    // class pointers
    System* System_pt;
    InfiniteSwitchSimulatedTempering* Isst_pt;
    Grid* Grid_pt;
    NptGrid* NptGrid_pt;
    SimulatedTempering* St_pt;

    // rotation steps
    void A_rot(Particle& particle, const double& h);
    void B_rot(Particle& particle, const double& h);
    void O_rot(Particle& particle);

    // NPT helper steps
    void A_NPT_2_box(const double& h);
    void B_NPT_box(const double& h);
    void O_NPT_box();
};

#endif
