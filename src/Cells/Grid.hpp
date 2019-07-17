#ifndef GRID_HPP
#define GRID_HPP

#include <omp.h>
#include <fstream>
#include <math.h>
#include <sstream>

#include "Generator.hpp"
#include "AverageObservable.hpp"
#include "Cells.hpp"
#include "Molecules.hpp"
#include "System.hpp"
#include "Array.hpp"
#include "SystemTemperature.hpp"
#include "RadialDistObservable.hpp"
#include "OrderObservable.hpp"
#include "SphereOrderObservable.hpp"
#include "SystemVolume.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                                GRID CLASS
// ------------------------------------------------------------------------- //
class Grid
{
public:
    // constructor
    Grid(const Matrix& Szero, const double& cut_off, Molecule* molecule_pt,
        const int& recf, const int& rect);
    // destructor
    ~Grid();
    // // function which calculates and returns the instant temperature
    // double get_instant_temperature();
    // // returns the currently held temperature
    // double get_temperature();
    // // updates the tracking objects by adding the current temperature to
    // its average
    void update_temperature();
    // fuction which updates the partice forces. It automatically checks if the
    // grid has changed sufficently to have to rebuild the grid. It then contineous
    // updateing the forces of all the particles using the linked lists from above,
    // looping over all the cells in the grid.
    void update_particle_forces(System* system_pt, Molecule* molecule_pt);
    // updates the positions from files given
    void add_file_initial_condition(Molecule* molecule_pt,
                                    const char* initial_pos_filename,
                                    const char* initial_mom_filename,
                                    const char* initial_box_filename);
    // set to calculate radial distribution
    void set_to_calculate_radial_dist(const double& rmin, const double& rmax,
                                      const int& N);
    void set_to_calculate_order_param(Molecule* molecule_pt,
                                      const double& part_rad);
    void set_to_calculate_sphere_order_param(Molecule* molecule_pt,
                                             const int& lmax,
                                             const double& cutoff);
    // access to observables
    SystemTemperature* Temperature_pt;
    SystemVolume* Volume_pt;

    RadialDistObservable* Radial_pt;
    OrderObservable* Order_pt;
    SphereOrderObservable* SphereOrder_pt;

    // access for grid dimensions
    Matrix S;
    Matrix Sp;

    // holer for exploding experiment
    bool break_experiment;

protected:
    // holder of the cut cut cut_off
    double cut_off_sq;
    // number of cells in each direction
    int number_of_cells_x;
    int number_of_cells_y;
    int number_of_cells_z;
    // number of particles in the grid
    unsigned number_of_particles;
    // bool check if we calculate radial dist
    bool with_radial_dist;
    // bool check if with order parameter
    bool with_order_param;
    bool with_sphere_order_param;
    // frequency to record the state at
    unsigned RecFreq;
    unsigned RecThres;

    // clear the forces and the potential of the molecule
    void clear_particle_forces(Molecule* molecule_pt);
    // function which maps a matrix index to a cell
    Cell* get_cell(const int& i, const int& j, const int& k);
    // check if the current grid satisfies the nessecary conditions, if it does
    // not then update the number of cells and return true boolean
    bool rebuild_grid_check();
    // rebuild a periodic grid
    void rebuild_periodic_grid(Molecule* molecule_pt);
    // update the position of all the particles on the grid
    void update_particles_on_grid();
    // returns the distance square between two particles
    Vector get_separation(const Particle& current_particle,
                          const Particle& neighbour_particle);
    // calculates and returns the momentum temperature
    double calculate_momentum_temp();
    // calulate and return the minimum image convention
    void calculate_min_image(Vector& r);
    // This function returns the non dimensional coordinate w.r.t the box
    // of a particle in the box.
    Vector get_box_coordinate(const Particle& particle);
    // reset the observables to record
    void reset_observables();

private:
    // controls the dimension of the grid
    // vector<double> L;
    // vector<double> Lp;
    // list of cells pointers
    vector<Cell*> cell_list;
    // tracking for average observables
    // AverageObservable* Temperature_pt;

    // add the particles from molecule object to the grid
    void add_particles_to_grid(Molecule* molecule_pt);
    // build a periodic grid
    void build_periodic_grid();
    // enforce the periodic boundary condition of the particle
    void enforce_periodic_particle_boundary_condition(Particle& particle);
    // delete and clear the grid
    void grid_clear();
    // returns the position in the grid of a particles based on
    // the particles position
    vector<int> get_cell_coordinate(Particle& particle_pt);
    // initialise the particles from the molecule on the grid
    void initialise_particles_on_grid(Molecule* molecule_pt);
    // set the neighbours such that they comunicate with each other
    // in the correct fashion. i.e so that it loops over the boundaries
    // of the matrix diagram
    void set_periodic_boundary_conditions();
    // initialise particles with random initial conditions, randomly uniform in
    // within the grid and with a gaussian distributed momentum.
    void set_random_particles_initial_condition(Molecule* molecule_pt);
    // computes the force in and respects the cut off radius. This is
    // a helper function because the iteration in the same cell is different
    // as we should not double count some distances
    void compute_force(System* system_pt, Molecule* molecule_pt,
                       Particle* current_particle, Particle* neighbour_particle);
    // read the initial position of the box from file
    void read_box_initial(const char* initial_box_filename);
    // read in the particle positions
    void read_position_initial(Molecule* molecule_pt,
                               const char* initial_pos_filename);
    // read in the particle momentums
    void read_momentum_initial(Molecule* molecule_pt,
                               const char* initial_mom_filename);
};

// ------------------------------------------------------------------------- //
//                                MOD HELPER
// ------------------------------------------------------------------------- //
static int mod(const int& a, const int& base)
{
    if(base != 0)
        return ((a % base) + base) % base;
    else
        return 0;
}

#endif
