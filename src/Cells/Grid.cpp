#include "Grid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                                GRID CLASS
// ------------------------------------------------------------------------- //
Grid::Grid(const Matrix& Szero, const double& cut_off,
           Molecule* molecule_pt) : S(2, 2), Sp(2, 2), break_experiment(false)
{
    if(Szero.size()[0] > 2)
    {
        S.resize(3, 3);
        Sp.resize(3, 3);
        number_of_cells_z = (int)floor(Szero(2, 2) / (2.0 * cut_off));

        if(number_of_cells_z == 0)
            number_of_cells_z = 1;

        number_of_neighbours = 10;
    }
    else
    {
        number_of_cells_z = 0;
        // set the number of neightbours for grids
        number_of_neighbours = 5;
    }

    // initialse the box dimensions
    S = Szero;

    // set the square of the cut off
    cut_off_sq = cut_off * cut_off;

    // calculate the number of boxes in each direction
    number_of_cells_x = (int)floor(S(0,0) / (2.0 * cut_off));
    number_of_cells_y = (int)floor(S(1,1) / (2.0 * cut_off));

    if(number_of_cells_x == 0)
        number_of_cells_x = 1;

    if(number_of_cells_y == 0)
        number_of_cells_y = 1;

    // build a periodic grid
    build_periodic_grid();

    // add all the particles to the grid
    initialise_particles_on_grid(molecule_pt);

    // set up tracking object to track the average
    // of the temperature
    Temperature_pt = new AverageObservable();

    // record the number of cells
    if(number_of_cells_z != 0)
        printf("Inital grid set as %d x %d x %d\n", number_of_cells_x,
                                                    number_of_cells_y,
                                                    number_of_cells_z);
    else
        printf("Inital grid set as %d x %d\n", number_of_cells_x,
                                               number_of_cells_y);
}

// destructor
Grid::~Grid()
{
    // delete the cells from the grid
    grid_clear();

    // delete all the average observables
    delete Temperature_pt;
}

// updates the positions from files given
void Grid::add_file_initial_condition(Molecule* molecule_pt,
                                      const char* initial_pos_filename,
                                      const char* initial_mom_filename,
                                      const char* initial_box_filename)
{
    // read in the particle positions
    read_position_initial(molecule_pt, initial_pos_filename);

    // read in the particle positions
    //read_momentum_initial(molecule_pt, initial_mom_filename);

    // read the box coordinates
    read_box_initial(initial_box_filename);

    // update the position of all particles
    update_particles_on_grid();
}

// add the particles from molecule object to the grid
void Grid::add_particles_to_grid(Molecule* molecule_pt)
{
    // loop over all the particles
    for(unsigned k=0; k<number_of_particles; k++)
    {
        // get the particle
        Particle * particle = &molecule_pt->particle(k);

        // enforce the periodc boundary condition on the particle
        enforce_periodic_particle_boundary_condition(*particle);

        // get the box coordinate based on particle positon
        vector<int> cell_coordinate = get_cell_coordinate(*particle);

        // assign this particle to the correct cell
        get_cell(cell_coordinate[0], cell_coordinate[1],
                 cell_coordinate[2])->assign_particle(particle);
    }
}

// build a periodic grid
void Grid::build_periodic_grid()
{
    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // build a new cell list
    for(int k=0; k<zend; k++)
        for(int j=0; j<number_of_cells_y; j++)
            for(int i=0; i<number_of_cells_x; i++)
                cell_list.push_back(new Cell);

    // set the boundary condition of the grid to be periodic
    set_periodic_boundary_conditions();
}

// clear the forces and the potential of the molecule
void Grid::clear_particle_forces(Molecule* molecule_pt)
{
    // loop over all the particles
    for(unsigned k=0; k<number_of_particles; k++)
    {
        // zero particle forces
        molecule_pt->particle(k).f.zero();

        if(molecule_pt->particle(k).rigid_body())
            molecule_pt->particle(k).tau.zero();
    }

    // clear the potential
    molecule_pt->potential() = 0.0;
}

// enforce the periodic boundary condition of the particle
void Grid::enforce_periodic_particle_boundary_condition(Particle& particle)
{
    Vector q_tilde = get_box_coordinate(particle);

    // enforce z
    if(number_of_cells_z != 0)
    {
        if(q_tilde(2) < 0.0)
            if(q_tilde(2) < -0.25)
            {
                printf("ERROR: particle detected in z: %1.3f\n", q_tilde(2));
                break_experiment = true;
            }
            else
                q_tilde(2) += 1.0;
        else if(q_tilde(2) > 1.0)
            if(q_tilde(2) > 1.25)
            {
                printf("ERROR: particle detected in z: %1.3f\n", q_tilde(2));
                break_experiment = true;
            }
            else
                q_tilde(2) -= 1.0;
    }

    // enforce y
    if(q_tilde(1) < 0.0)
        if(q_tilde(1) < -0.25)
        {
            printf("ERROR: particle detected in y: %1.3f\n", q_tilde(1));
            break_experiment = true;
        }
        else
            q_tilde(1) += 1.0;
    else if(q_tilde(1) > 1.0)
        if(q_tilde(1) > 1.25)
        {
            printf("ERROR: particle detected in y: %1.3f\n", q_tilde(1));
            break_experiment = true;
        }
        else
            q_tilde(1) -= 1.0;

    // enforce x
    if(q_tilde(0) < 0.0)
        if(q_tilde(0) < -0.25)
        {
            printf("ERROR: particle detected in x: %1.3f\n", q_tilde(0));
            break_experiment = true;
        }
        else
            q_tilde(0) += 1.0;
    else if(q_tilde(0) > 1.0)
        if(q_tilde(0) > 1.25)
        {
            printf("ERROR: particle detected in x: %1.3f\n", q_tilde(0));
            break_experiment = true;
        }
        else
            q_tilde(0) -= 1.0;

    // recalculate position
    particle.q = S * q_tilde;
}

// delete and clear the grid
void Grid::grid_clear()
{
    // clear the current cell list
    cell_list.clear();
}
//
// // function which maps a matrix index to a cell
// Cell* Grid::get_cell(const int& i, const int& j)
// {
//     // define helpers for number of cells
//     // in each direction
//     int nx = number_of_cells_x;
//     int ny = number_of_cells_y;
//
//     return cell_list[mod(i, nx) + mod(j, ny) * nx];
// }
//

// function which maps a matrix index to a cell
Cell* Grid::get_cell(const int& i, const int& j, const int& k)
{
    // define helpers for number of cells
    // in each direction
    int nx = number_of_cells_x;
    int ny = number_of_cells_y;
    int nz = number_of_cells_z;

    return cell_list[mod(i, nx) + mod(j, ny) * nx + mod(k, nz) * nx * ny];
}

// returns the number of grid coordinates
unsigned Grid::get_ncoord()
{
    return 3;
}


// This function returns the non dimensional coordinate w.r.t the box
// of a particle in the box.
Vector Grid::get_box_coordinate(const Particle& particle)
{

  // rescale the coordiinate
  Vector q_tilde = S.inv() * particle.q;

  // return the coordinates
  return q_tilde;
}

// returns the position in the grid of a particles based on
// the particles position
std::vector<int> Grid::get_cell_coordinate(Particle& particle)
{
    // make a vector of coordinates
    std::vector<int> cell_coordinate(3, 0);
    Vector q_tilde = get_box_coordinate(particle);

    cell_coordinate[0] = (int)floor(q_tilde(0) * number_of_cells_x);
    cell_coordinate[1] = (int)floor(q_tilde(1) * number_of_cells_y);

    if(number_of_cells_z != 0)
        cell_coordinate[2] = (int)floor(q_tilde(2) * number_of_cells_z);

    // return the vector with the coordinates
    return cell_coordinate;
}

// returns the distance square between two particles
Vector Grid::get_separation(const Particle& current_particle,
                            const Particle& neighbour_particle)
{
   // Calculate separation
    Vector r = neighbour_particle.q - current_particle.q;

    // update difference so that follows mini image
    calculate_min_image(r);

    // return the separation
    return r;
}

// calulate and return the minimum image convention
void Grid::calculate_min_image(Vector& r)
{
    // enforce minimum image convention
    if(r(0) > 0.5 * S(0,0))
        r(0) -= S(0,0);
    else if(r(0) <= -0.5 * S(0,0))
        r(0) += S(0,0);

    // enforce minimum image convention
    if(r(1) > 0.5 * S(1,1))
        r(1) -= S(1,1);
    else if(r(1) <= -0.5 * S(1,1))
        r(1) += S(1,1);

    if(number_of_cells_z != 0)
    {
        if(r(2) > 0.5 * S(2, 2))
            r(2) -= S(2, 2);
        else if(r(2) <= -0.5 * S(2, 2))
            r(2) += S(2, 2);
    }
}

// calculates the translational shift of the box
// as a function of y which is the position of a
// particle
// double Grid::get_translational_shift_x(const Vector& r)
// {
//     if(number_of_cells_z != 0)
//         return r(1) * S(0, 1) / S(1, 1) + r(2) * S(0, 2) / S(2, 2);
//     else
//         return r(1) * S(0, 1) / S(1, 1);
// }
//
// double Grid::get_translational_shift_y(const Vector& r)
// {
//     if(number_of_cells_z != 0)
//         return r(2) * S(1, 2) / S(2, 2);
//     else
//         return 0.0;
// }

// initialise the particles from the molecule on the grid
void Grid::initialise_particles_on_grid(Molecule* molecule_pt)
{
    // set the number of particles
    number_of_particles = molecule_pt->nparticle();

    // set random initial conditions
    set_random_particles_initial_condition(molecule_pt);

    // add the particles to the grid
    add_particles_to_grid(molecule_pt);
}

// check if the current grid satisfies the nessecary conditions, if it does
// not then update the number of cells and return true boolean
bool Grid::rebuild_grid_check()
{
    // boolean if the grid needs to change
    bool need_to_change_grid = false;

    // calculate the correct number of cells
    int new_number_of_cells_x = (int)floor(S(0,0) / (2.0 * sqrt(cut_off_sq)));
    int new_number_of_cells_y = (int)floor(S(1,1) / (2.0 * sqrt(cut_off_sq)));
    int new_number_of_cells_z = 0;

    if(new_number_of_cells_x == 0)
        new_number_of_cells_x = 1;

    if(number_of_cells_y == 0)
        new_number_of_cells_y = 1;

    if(number_of_cells_z != 0)
    {
        new_number_of_cells_z = (int)floor(S(2,2) / (2.0 * sqrt(cut_off_sq)));

        if(new_number_of_cells_x == 0)
            new_number_of_cells_x = 1;

        // perform the check in the y direction
        if(new_number_of_cells_z != number_of_cells_z)
        {
            // check that cell is sensible
            if(new_number_of_cells_z < 0 ||
                new_number_of_cells_z > number_of_cells_z + 2)
                {
                    printf("\nERROR: cell is exploding new cell would be: %d x %d x %d\n\n",
                            new_number_of_cells_x, new_number_of_cells_y,
                            new_number_of_cells_z);

                    // set break experiment to true
                    break_experiment = true;
                    // exit(-1);
                }
            need_to_change_grid = true;
            number_of_cells_z = new_number_of_cells_z;
        }
    }

    // perform check in x direction
    if(new_number_of_cells_x != number_of_cells_x)
    {
        // check that cell is sensible
        if(new_number_of_cells_x < 0 ||
            new_number_of_cells_x > number_of_cells_x + 2)
            {
                printf("\nERROR: cell is exploding new cell would be: %d x %d x %d\n\n",
                        new_number_of_cells_x, new_number_of_cells_y,
                        new_number_of_cells_z);

                // set break experiment to true
                break_experiment = true;
                // exit(-1);
            }

        need_to_change_grid = true;
        number_of_cells_x = new_number_of_cells_x;
    }

    // perform the check in the y direction
    if(new_number_of_cells_y != number_of_cells_y)
    {
        // check that cell is sensible
        if(new_number_of_cells_y < 0 ||
            new_number_of_cells_y > number_of_cells_y + 2)
            {
                printf("\nERROR: cell is exploding new cell would be: %d x %d x %d\n\n",
                        new_number_of_cells_x, new_number_of_cells_y,
                        new_number_of_cells_z);
                // set break experiment to true
                break_experiment = true;
                // exit(-1);
            }
        need_to_change_grid = true;
        number_of_cells_y = new_number_of_cells_y;
    }

    // if(need_to_change_grid)
    //     printf("New Grid is: %d x %d x %d\n", number_of_cells_x,
    //             number_of_cells_y, number_of_cells_z);

    // return the result of the check
    return need_to_change_grid;
}

// rebuild a periodic grid
void Grid::rebuild_periodic_grid(Molecule* molecule_pt)
{
    // delete the current grid
    grid_clear();

    // build a periodic grid again
    build_periodic_grid();

    // add the particles to the grid
    add_particles_to_grid(molecule_pt);
}

// set the neighbours such that they comunicate with each other
// in the correct fashion. i.e so that it loops over the boundaries
// of the matrix diagram
void Grid::set_periodic_boundary_conditions()
{
    // There are 5 neighbours in total and should be allocated as:
    // +---+--------+----------+
    // |   |        | i+1, j-1 |
    // +---+--------+----------+
    // |   | i, j   | i+1, j   |
    // +---+--------+----------+
    // |   | i, j+1 | i+1, j+1 |
    // +---+--------+----------+
    for(int j=0; j<number_of_cells_y; j++)
        for(int i=0; i<number_of_cells_x; i++)
        {
            if(number_of_cells_z != 0)
            {
                for(int k=0; k<number_of_cells_z; k++)
                {
                    // set the neightbours
                    get_cell(i, j, k)->
                                set_neighbour(0, get_cell(i+1, j-1, k));
                    get_cell(i, j, k)->
                                set_neighbour(1, get_cell(i, j, k));
                    get_cell(i, j, k)->
                                set_neighbour(2, get_cell(i+1, j, k));
                    get_cell(i, j, k)->
                                set_neighbour(3, get_cell(i, j+1, k));
                    get_cell(i, j, k)->
                                set_neighbour(4, get_cell(i+1, j+1, k));
                    get_cell(i, j, k)->
                                set_neighbour(5, get_cell(i+1, j-1, k+1));
                    get_cell(i, j, k)->
                                set_neighbour(6, get_cell(i, j, k+1));
                    get_cell(i, j, k)->
                                set_neighbour(7, get_cell(i+1, j, k+1));
                    get_cell(i, j, k)->
                                set_neighbour(8, get_cell(i, j+1, k+1));
                    get_cell(i, j, k)->
                                set_neighbour(9, get_cell(i+1, j+1, k+1));
                }
            }
            else
            {
                // set the neightbours
                get_cell(i, j, 0)->
                            set_neighbour(0, get_cell(i+1, j-1, 0));
                get_cell(i, j, 0)->
                            set_neighbour(1, get_cell(i, j, 0));
                get_cell(i, j, 0)->
                            set_neighbour(2, get_cell(i+1, j, 0));
                get_cell(i, j, 0)->
                            set_neighbour(3, get_cell(i, j+1, 0));
                get_cell(i, j, 0)->
                            set_neighbour(4, get_cell(i+1, j+1, 0));
            }
        }
}

// initialise particles with random initial conditions, randomly uniform in
// within the grid and with a gaussian distributed momentum.
void Grid::set_random_particles_initial_condition(Molecule* molecule_pt)
{
  // derefernce the number of particles
  double beta = 1.0 / molecule_pt->kt();

  // set up random normal generator
  NormalGenerator normal_momentum(0.0, 1.0, 4563);

  int particles_in_dir =
        (int)ceil(pow(number_of_particles,
                      1.0 / (2 + (number_of_cells_z != 0))));

  // hexagonal separation
  double separation_in_x = S(0,0) / double(particles_in_dir + 1.0);
  double separation_in_y = S(1,1) / double(particles_in_dir + 1.0);
  double separation_in_z = 0.0;

  if(number_of_cells_z != 0)
    separation_in_z = S(2, 2) / double(particles_in_dir + 1.0);

  // particle counter
  unsigned pc = 0;

  int zend = 1;
    if(number_of_cells_z != 0)
        zend = particles_in_dir;

  // loop over each direction
  for(int k=0; k<zend; k++)
    for(int j=0; j<particles_in_dir; j++)
        for(int i=0; i<particles_in_dir; i++)
        {
            if(pc < number_of_particles)
            {
                // derefernce the particle
                Particle* particle = &molecule_pt->particle(pc);

                particle->q(0) = double(i+1) * separation_in_x;
                particle->q(1) = double(j+1) * separation_in_y;

                if(number_of_cells_z != 0)
                    particle->q(2) = double(k+1) * separation_in_z;

                double mx = particle->m(0,0);
                double my = particle->m(1,1);
                double mz = 0.0;

                if(number_of_cells_z != 0)
                    mz = particle->m(2,2);

                particle->p(0) = sqrt(mx / beta) * normal_momentum();
                particle->p(1) = sqrt(my / beta) * normal_momentum();

                if(number_of_cells_z != 0)
                    particle->p(2) = sqrt(mz / beta) * normal_momentum();

                // if particle rotates set all the rotations
                if(particle->rigid_body() != false)
                {
                    double alpha = 0.0;

                    // rotate only uneven particles
                    if(j % 2 == 0)
                    {
                        if(i % 2 != 0 && i != 0)
                            alpha = M_PI/3.0;
                    }
                    else if(j % 2 == 1)
                    {
                        if(i % 2 == 0)
                            alpha = M_PI/3.0;
                    }

                    // set rotation
                    RotMatrix Rot(alpha);
                    particle->Q = Rot;
                }
            }

            // increment the particle counter
            pc++;
        }
}

// update the position of all the particles on the grid
void Grid::update_particles_on_grid()
{
    // make parameters
    ListNode* cell_conductor = NULL;
    ListNode* cell_next = NULL;

    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // loop over all the cells
    for(int k=0; k<zend; k++)
        for(int j=0; j<number_of_cells_y; j++)
            for(int i=0; i<number_of_cells_x; i++)
            {
                // get the conductor for this cell
                cell_conductor = get_cell(i, j, k)->get_particle_list_head();

                // loop over all the particles in this cell
                while(cell_conductor != NULL)
                {
                    // derefernce the next particle
                    // (need to hold this in case we move particle such that
                    // the list is not changed)
                    cell_next = cell_conductor->next;

                    // dereference the partilce
                    Particle* particle = cell_conductor->particle;

                    // if we are in boundary cell then enforce the
                    // periodic boundary condition
                    if(i == 0 or i == number_of_cells_y-1)
                        enforce_periodic_particle_boundary_condition(*particle);
                    else if(j == 0 or j == number_of_cells_x-1)
                        enforce_periodic_particle_boundary_condition(*particle);
                    else if(number_of_cells_z != 0 and k == 0 or k == number_of_cells_z-1)
                        enforce_periodic_particle_boundary_condition(*particle);

                    // get the potentially new cell cordinates
                    vector<int> cell_coordinate =
                                            get_cell_coordinate(*particle);

                    //printf("Cell coordinate %d %d %d\n", cell_coordinate[0],cell_coordinate[1],cell_coordinate[2]);

                    // check if particle has moved cell
                    if(cell_coordinate[0] != i
                        or cell_coordinate[1] != j or cell_coordinate[2] != k)
                    {
                        // dereference for readability
                        int new_i = cell_coordinate[0];
                        int new_j = cell_coordinate[1];
                        int new_k = cell_coordinate[2];

                        // move the particle to the new cell
                        get_cell(i, j, k)->
                                move_list_node(cell_conductor,
                                               get_cell(new_i, new_j, new_k));
                    }

                    // step the conductor forward
                    cell_conductor = cell_next;
                }
            }
}

// function which calculates and returns the
// volume of the current cell
double Grid::get_instant_temperature()
{
  return Temperature_pt->get_instant();
}

// returns the currently held temperature
double Grid::get_temperature()
{
  return Temperature_pt->get_average();
}

// updates the tracking objects by adding the
// current temperature to its average
void Grid::update_temperature()
{
  // calculate the current temperature
  double temperature = calculate_momentum_temp();

  // add this value to the tracking average
  Temperature_pt->observe(temperature);
}

// calculates and returns the momentum temperature
double Grid::calculate_momentum_temp()
{
    // initialise the temperature
    double temperature = 0.0;
    double rot_temperature = 0.0;

    // helpers
    Particle* particle = NULL;
    ListNode* cell_conductor = NULL;

    double N = 1.0;

    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // loop over all the cells
    for(int k=0; k<zend; k++)
        for(int j=0; j<number_of_cells_y; j++)
            for(int i=0; i<number_of_cells_x; i++)
            {
                // get the conductor for this cell
                cell_conductor = get_cell(i, j, k)->get_particle_list_head();

                // loop over all the particles in this cell
                while(cell_conductor != NULL)
                {
                    // dereference the partilce
                    particle = cell_conductor->particle;

                    temperature = particle->p.dot(particle->m * particle->p) / N
                                + (N - 1.0) / N * temperature;

                    //if particle rotates add to rotational temperature
                    if(particle->rigid_body())
                    {
                        double pi = particle->pi(0,0);
                        double I =  particle->I(0,0);

                        rot_temperature = pi * pi / (I * N)
                                    + (N-1.0)/N * rot_temperature;
                    }

                    // increment the number of particles
                    N++;

                    // step the conductor forward
                    cell_conductor = cell_conductor->next;
                }
            }

    // total temperature
    double total_temp = 0.0;

    // multiply by the correct constant
    if(rot_temperature != 0.0)
        total_temp = 1.0 / 3.0 * (temperature + rot_temperature);
    else
        total_temp = 1.0 / (2.0 + (number_of_cells_z != 0)) * temperature;

    // return the temperature
    return total_temp;
}

// computes the force in and respects the cut off radius. This is
// a helper function because the iteration in the same cell is different
// as we should not double count some distances
void Grid::compute_force(System* system_pt, Molecule* molecule_pt,
                         Particle* current_particle,
                         Particle* neighbour_particle)
{
    // calculate the square of the distance between particles
    Vector r = get_separation(*current_particle, *neighbour_particle);

    double rsq = r.l22();

    // chcek the cutoff criterion
    if(rsq < cut_off_sq)
        // Calculate force and potential
        system_pt->compute_force(molecule_pt, current_particle,
                                 neighbour_particle, std::sqrt(rsq), r);
}

// fuction which updates the partice forces. It updates the forces of all
// the particles using the linked lists from above, looping over all the
// cells in the grid.
void Grid::update_particle_forces(System* system_pt, Molecule* molecule_pt)
{
    // make sure taht the particles are
    // on the grid and that the forces are set to zero
    #pragma omp single
    {
        // update the particles on the grid
        update_particles_on_grid();

        // clear all the particle forces
        clear_particle_forces(molecule_pt);
    }

    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

  // loop over all the boxes to calculate the forces
#pragma omp parallel for default(shared) schedule(dynamic) collapse(3) \
      firstprivate(zend, number_of_neighbours, system_pt)
    for(int k=0; k<zend; k++)
    {
        for(int j=0; j<number_of_cells_y; j++)
        {
            for(int i=0; i<number_of_cells_x; i++)
            {
                for(unsigned n=0; n<number_of_neighbours; n++)
                {
                    // get the current cell
                    Cell* current_cell = get_cell(i, j, k);

                    // get the current cells conductor
                    ListNode* current_conductor =
                                    current_cell->get_particle_list_head();

                    // get the neighbour cell
                    Cell* neighbour_cell = current_cell->get_neighbour(n);

                    // reset newton iterator
                    unsigned current_newton_iterator = 0;

                    // loop over particles in current cell
                    while(current_conductor != NULL)
                    {
                        // conductor over neighbour cell
                        ListNode* neighbour_conductor =
                                    neighbour_cell->get_particle_list_head();

                        // get the particle pointer
                        Particle* current_particle =
                                    current_conductor->particle;

                        // reset newton iterator
                        unsigned neighbour_newton_iterator = 0;

                        // loop over particles in neighbour cell
                        while(neighbour_conductor != NULL)
                        {
                            // apply Newton iterators if we are in the same cell
                            if(current_cell == neighbour_cell)
                            {
                                if(current_newton_iterator <=
                                                neighbour_newton_iterator)
                                {
                                    // get the neighbour particle
                                    Particle* neighbour_particle =
                                                neighbour_conductor->particle;

                                    // don't compute the force for the same
                                    // particles
                                    if(current_particle != neighbour_particle)
                                    {
                                        compute_force(system_pt, molecule_pt,
                                                      current_particle,
                                                      neighbour_particle);
                                    }
                                }

                                // iterate the newton current newton iterator
                                neighbour_newton_iterator++;
                            }
                            else
                            {
                                // get the neighbour particle
                                Particle* neighbour_particle =
                                                neighbour_conductor->particle;

                                // compute the force
                                compute_force(system_pt, molecule_pt,
                                              current_particle,
                                              neighbour_particle);
                            }

                            // itreate the current conductor
                            neighbour_conductor = neighbour_conductor->next;

                        } // end of loop over neighbour_conductor

                        // iterate the newton current newton iterator
                        current_newton_iterator++;

                        // itreate the current conductor
                        current_conductor = current_conductor->next;

                    } // end of loop over current_conductor
                } // end of loop over n, number of neighbours
            } // end of loop over i, x-direction
        } // end of loop over j, y-direction
    } // end of loop over k, z-direction
}

 // read the initial position of the box from file
 void Grid::read_box_initial(const char* initial_box_filename)
 {
     // loop over the box file stream
     std::ifstream input(initial_box_filename);
     std::string line;

    // loop over all the lines in the box
    while(std::getline(input, line))
    {
        std::stringstream lineStream(line);
        std::string line_stuff;
        std::vector<double> parsed_row;

        while(std::getline(lineStream, line_stuff, ','))
        {
            double test = strtod(line_stuff.c_str(), NULL);
            parsed_row.push_back(test);
        }

       for(unsigned i=0; i<parsed_row.size();i++)
       {
           if(i<3)
           {
               if(i == 0)
                    S(0,0) = parsed_row[i];
                else if(i == 1)
                    S(0,1) = parsed_row[i];
                else
                    S(1,1) = parsed_row[i];
           }
           else
           {
               if(i == 0)
                    Sp(0,0) = parsed_row[i];
                else if(i == 1)
                    S(0,1) = parsed_row[i];
                else
                    S(1,1) = parsed_row[i];
           }
       }
    }
}

// read in the particle positions
void Grid::read_position_initial(Molecule* molecule_pt,
                                 const char* initial_pos_filename)
{
   // loop over the box file stream
   std::ifstream input(initial_pos_filename);
   std::string line;

   Particle* particle = NULL;
   unsigned k = 0;

  // loop over all the lines in the box
  while(std::getline(input, line))
  {
      std::stringstream lineStream(line);
      std::string line_stuff;
      std::vector<double> parsed_row;

      while(std::getline(lineStream, line_stuff, ','))
      {
          double test = strtod(line_stuff.c_str(), NULL);
          parsed_row.push_back(test);
      }

      particle = &molecule_pt->particle(k);

      // set position
      particle->q(0) = parsed_row[0];
      particle->q(1) = parsed_row[1];

      // if we rotates
      if(particle->rigid_body())
      {
          particle->Q(0, 0) = parsed_row[2];
          particle->Q(0, 1) = parsed_row[3];
          particle->Q(1, 0) = parsed_row[4];
          particle->Q(1, 1) = parsed_row[5];
      }

      // increment particle counter
      k++;

     }
 }

// read in the particle momentums
void Grid::read_momentum_initial(Molecule* molecule_pt,
                                 const char* initial_mom_filename)
{
    // loop over the box file stream
    std::ifstream input(initial_mom_filename);
    std::string line;

    Particle* particle = NULL;
    unsigned k = 0;

   // loop over all the lines in the box
   while(std::getline(input, line))
   {
       std::stringstream lineStream(line);
       std::string line_stuff;
       std::vector<double> parsed_row;

       while(std::getline(lineStream, line_stuff, ','))
       {
           double test = strtod(line_stuff.c_str(), NULL);
           parsed_row.push_back(test);
       }

       particle = &molecule_pt->particle(k);

       // set momentum
       particle->p(0) = parsed_row[0];
       particle->p(1) = parsed_row[1];

       // if we rotates
       if(particle->rigid_body())
       {
           particle->pi(0,0) = parsed_row[2];
       }

       // increment particle counter
       k++;
      }
}
