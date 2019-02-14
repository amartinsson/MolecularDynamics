#ifndef GRID_HPP
#define GRID_HPP

#include <omp.h>
#include <fstream>
#include <sstream>

#include "Generator.hpp"
#include "AverageObservable.hpp"
#include "Cells.hpp"
#include "Molecules.hpp"
#include "System.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                                GRID CLASS
// ------------------------------------------------------------------------- //
class Grid
{
public:
    // constructor
    Grid(const double& a_x, const double& b_x, const double& b_y,
         const double& cut_off, Molecule* molecule_pt);
    // destructor
    ~Grid();
    // access for grid dimensions
    double* L_pt(const unsigned& i);
    // access for grid dimension momentum
    double* Lp_pt(const unsigned& i);
    // function which calculates and returns the instant temperature
    double get_instant_temperature();
    // returns the number of grid coordinates
    unsigned get_ncoord();
    // returns the currently held temperature
    double get_temperature();
    // updates the tracking objects by adding the current temperature to
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
protected:
    // holder of the cut cut cut_off
    double cut_off_sq;
    // number of cells in each direction
    int number_of_cells_x;
    int number_of_cells_y;
    // number of neightbours for each cell
    unsigned number_of_neighbours;
    // number of particles in the grid
    unsigned number_of_particles;

    // clear the forces and the potential of the molecule
    void clear_particle_forces(Molecule* molecule_pt);
    // function which maps a matrix index to a cell
    Cell* get_cell(const int& i, const int& j);
    // check if the current grid satisfies the nessecary conditions, if it does
    // not then update the number of cells and return true boolean
    bool rebuild_grid_check();
    // rebuild a periodic grid
    void rebuild_periodic_grid(Molecule* molecule_pt);
    // update the position of all the particles on the grid
    void update_particles_on_grid();
    // returns the distance square between two particles
    vector<double> get_distance_square(Particle* current_particle,
                                       Particle* neighbour_particle);
    // calculates and returns the momentum temperature
    double calculate_momentum_temp();

private:
    // controls the dimension of the grid
    vector<double> L;
    vector<double> Lp;
    // list of cells pointers
    vector<Cell*> cell_list;
    // tracking for average observables
    AverageObservable* Temperature_pt;

    // add the particles from molecule object to the grid
    void add_particles_to_grid(Molecule* molecule_pt);
    // build a periodic grid
    void build_periodic_grid();
    // enforce the periodic boundary condition of the particle
    void enforce_periodic_particle_boundary_condition(Particle* part_pt);
    // delete and clear the grid
    void grid_clear();
    // returns the position in the grid of a particles based on
    // the particles position
    vector<int> get_cell_coordinate(Particle* particle_pt);
    // calculates the translational shift of the box
    // as a function of y which is the position of a
    // particle
    double get_translational_shift(const double& y);
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
    return ((a % base) + base) % base;
}

#endif
