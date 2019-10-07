#ifndef NPTGRID_HPP
#define NPTGRID_HPP

#include "Grid.hpp"
#include "SystemPressure.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                               NPT GRID CLASS
// ------------------------------------------------------------------------- //
class NptGrid : public Grid
{
public:
    NptGrid(const Matrix& Szero,
            const double& cut_off, Molecule* molecule_pt,
            const double& mass, const double& target_press,
            const int& recf, const int& rect);
    // destructor
    ~NptGrid();
    // update both the pressure and the temperature variables -
    // this needs to be called at the end of every integration step
    // and will update the pressure and temperature
    void update_pressure_temperature();
    // This function takes in a 2 dimensional vector and returns the
    // box normalised vector r_tilde
    // Vector get_box_min_image_sep(const Particle& current_particle,
    //                              const Particle& neighbour_particle);
    // function which sets the position of a particle to the invariant
    // position given by q tilde
    void set_box_coordinate(const Vector& q_tilde, Particle& particle);
    // This function calculates the non dimensional momentum coordinate w.r.t
    // to the box of the particle given
    Vector get_box_momentum(const Particle& particle);
    // function which sets the position of a particle to the invariant
    // position given by q tilde
    void set_box_momentum(const Vector& p_tilde, Particle& particle);
    // get the mass of the grid
    double get_mass();
    double get_target_pressure() const {return target_pressure;}

    void reset_target_pressure(const double& new_pressure,
    const int& recf, const int& rect) {

        // set the new target pressure
        target_pressure = new_pressure;

        // reset the observables
        RecThres = rect;
        RecFreq = recf;

        Volume_pt->reset(RecFreq, RecThres);
        Temperature_pt->reset(RecFreq, RecThres);
        Pressure_pt->reset(RecFreq, RecThres);

        Radial_pt->reset(RecFreq, RecThres);
    };
    // enforces the relaative positoon of the particles
    void enforce_constant_relative_particle_pos(const Matrix& Sold);
    // returns the box force in direction i
    void update_kinetic_gradient();
    // returns the virial for the ith direction
    // double get_virial(const unsigned& i){return box_grad_potential[i];}
    // fuction which updates the particle forces. It automatically checks if the
    // grid has changed sufficently to have to rebuild the grid. It then contineous
    // updateing the forces of all the particles using the linked lists from above,
    // looping over all the cells in the grid. It also updates the pressure and
    // the temperature by calculating the dot product of q and f for all the
    // particles
    void update_particle_forces(System* system_pt, Molecule* molecule_pt);

    // Stores the virial function
    Matrix virial;
    Matrix nablaK;
    // tracking for average observables
    SystemPressure* Pressure_pt;

private:

    // target pressure
    double target_pressure;
    // holder for the mass of the grid --
    // this is needed for NPT integration
    double box_mass;
    // holder for box force
    // std::vector<double> box_grad_potential;
    double box_grad_00;
    double box_grad_01;
    double box_grad_11;

    double box_grad_02;
    double box_grad_12;
    double box_grad_22;

    //double* box_grad_potential;
    Vector momentum_sq;
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
    // void update_pressure();
    // function which updates the volume
    // void update_volume();
    // enforce the constraint that the relative distances
    // of the particle position and momentum cannot change.
    void enforce_relative_particle(const Matrix& Sold);
    // reset the observables for the next step
    // void reset_observables();

};

#endif
