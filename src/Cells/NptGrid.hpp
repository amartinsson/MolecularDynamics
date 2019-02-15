#ifndef NPTGRID_HPP
#define NPTGRID_HPP

#include "Grid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                               NPT GRID CLASS
// ------------------------------------------------------------------------- //
class NptGrid : public Grid
{
public:
    NptGrid(const double& a_x, const double& b_x, const double& b_y,
            const double& cut_off, Molecule* molecule_pt,
            const double& mass, const double& target_press);
    // destructor
    ~NptGrid();
    // update both the pressure and the temperature variables -
    // this needs to be called at the end of every integration step
    // and will update the pressure and temperature
    void update_pressure_temperature();
    // This function returns the non dimensional coordinate w.r.t the box
    // of a particle in the box.
    vector<double> get_box_coordinate(Particle* particle_pt);
    // This function takes in a 2 dimensional vector and returns the
    // box normalised vector r_tilde
    vector<double> get_box_min_image_sep(Particle* current_particle,
                                         Particle* neighbour_particle);
    // function which sets the position of a particle to the invariant
    // position given by q tilde
    void set_box_coordinate(vector<double>& q_tilde, Particle* particle_pt);
    // This function calculates the non dimensional momentum coordinate w.r.t
    // to the box of the particle given
    vector<double> get_box_momentum(Particle* particle_pt);
    // function which sets the position of a particle to the invariant
    // position given by q tilde
    void set_box_momentum(std::vector<double>& p_tilde, Particle* particle_pt);
    // get the mass of the grid
    double get_mass();
    // function which calculates and returns the
    // volume of the current cell
    double get_instant_volume();
    // function which returns the instantaneous pressure
    double get_instant_pressure();
    // function which returns the cahced volume  which is calculated when
    // temperature and pressure are updated
    double get_volume();
    // returns the currently held pressure
    double get_pressure();
    // enforces the relaative positoon of the particles
    void enforce_constant_relative_particle_pos(const vector<double>& L_old);
    // returns the box force in direction i
    double get_accumulated_momentum(const unsigned& i){return momentum_sq[i];}
    // updates the accumulated momentum
    void update_accumulted_momentum();
    // returns the virial for the ith direction
    double get_virial(const unsigned& i){return box_grad_potential[i];}
    // fuction which updates the particle forces. It automatically checks if the
    // grid has changed sufficently to have to rebuild the grid. It then contineous
    // updateing the forces of all the particles using the linked lists from above,
    // looping over all the cells in the grid. It also updates the pressure and
    // the temperature by calculating the dot product of q and f for all the
    // particles
    void update_particle_forces(System* system_pt, Molecule* molecule_pt);

private:
    // tracking for average observables
    AverageObservable* Pressure_pt;
    AverageObservable* Volume_pt;
    // target pressure
    double target_pressure;
    // holder for the mass of the grid --
    // this is needed for NPT integration
    double box_mass;
    // holder for box force
    std::vector<double> box_grad_potential;
    double box_grad_zero;
    double box_grad_one;
    double box_grad_two;
    //double* box_grad_potential;
    std::vector<double> momentum_sq;
    // reducto helpers for caluclating the force
    double red_help_grad_ax;
    double red_help_grad_bx;
    double red_help_grad_by;

    // computes the force in and respects the cut off radius. This is
    // a helper function because the iteration in the same cell is different
    // as we should not double count some distances
    void compute_force(System* system_pt, Molecule* molecule_pt,
                       Particle* current_particle,
                       Particle* neighbour_particle);
    // function which updates the pressure
    void update_pressure();
    // function which updates the volume
    void update_volume();
    // enforce the constraint that the relative distances
    // of the particle position and momentum cannot change.
    void enforce_relative_particle(const vector<double> L_old);
    // function which calculates the pressure
    double calculate_pressure();
};

#endif
