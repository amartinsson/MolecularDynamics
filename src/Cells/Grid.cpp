#include "Grid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                                GRID CLASS
// ------------------------------------------------------------------------- //
Grid::Grid(const double& a_x, const double& b_x, const double& b_y,
           const double& cut_off, Molecule* molecule_pt)
{
    // initialise the box vector;
    L.resize(3);
    Lp.resize(3);

    // initialise the box dimensions
    L[0] = a_x;
    L[1] = b_x;
    L[2] = b_y;

    // initialise the box momentum to zero
    for(unsigned i=0; i<3; i++)
    Lp[i] = 0.0;

    // set the square of the cut off
    cut_off_sq = cut_off * cut_off;

    // calculate the number of boxes in each direction
    number_of_cells_x = (int)floor(L[0] / (2.0 * cut_off));
    number_of_cells_y = (int)floor(L[2] / (2.0 * cut_off));

    // build a periodic grid
    build_periodic_grid();

    // set the number of neightbours for grids
    number_of_neighbours = 5;

    // add all the particles to the grid
    initialise_particles_on_grid(molecule_pt);

    // set up tracking object to track the average
    // of the temperature
    Temperature_pt = new AverageObservable();
}

// destructor
Grid::~Grid()
{
    // delete the cells from the grid
    grid_clear();

    // clear the memory of the vectors
    L.clear();
    Lp.clear();

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
    // derefernce the number of particles
    std::vector<int> cell_coordinate;
    Particle* particle_k = NULL;

    // loop over all the particles
    for(unsigned k=0; k<number_of_particles; k++)
    {
        // get the particle
        particle_k = molecule_pt->particle_pt(k);

        // enforce the periodc boundary condition on the particle
        enforce_periodic_particle_boundary_condition(particle_k);

        // get the box coordinate based on particle positon
        cell_coordinate = get_cell_coordinate(particle_k);

        // assign this particle to the correct cell
        get_cell(cell_coordinate[0], cell_coordinate[1])->
                                                assign_particle(particle_k);
    }
}

// build a periodic grid
void Grid::build_periodic_grid()
{
    // build a new cell list
    for(int j=0; j<number_of_cells_y; j++)
    for(int i=0; i<number_of_cells_x; i++)
    cell_list.push_back(new Cell);

    // set the boundary condition of the grid to be periodic
    set_periodic_boundary_conditions();
}

// clear the forces and the potential of the molecule
void Grid::clear_particle_forces(Molecule* molecule_pt)
{
    // dereference the number of particles
    unsigned dim = molecule_pt->dim();

    // make a holder for particle k
    Particle* particle_k = NULL;

    // loop over all the particles
    for(unsigned k=0; k<number_of_particles; k++)
    {
        // dereference tha particle
        particle_k = molecule_pt->particle_pt(k);

        // loop over the dimensions
        for(unsigned j=0; j<dim; j++)
        {
            *particle_k->f_pt(j) = 0.0;
        }

        // if the particle rotates - reset torque
        if(particle_k->rigid_body() != false)
        *particle_k->tau_pt(0) = 0.0;
    }

    // clear the potential
    *molecule_pt->potential_pt() = 0.0;
}

// enforce the periodic boundary condition of the particle
void Grid::enforce_periodic_particle_boundary_condition(Particle* part_pt)
{
    // dereference the position
    double x = *part_pt->q_pt(0);
    double y = *part_pt->q_pt(1);

    // enforce periodicity in y direction
    if(y < 0.0)
    y += L[2];
    else if(y > L[2])
    y -= L[2];

    // get the translational shift of this particle
    double trans_shift = get_translational_shift(y);

    // enforce periodicity in x direction
    if(x < trans_shift)
    x += L[0];
    else if (x > L[0] + trans_shift)
    x -= L[0];

    // set the new values
    *part_pt->q_pt(0) = x;
    *part_pt->q_pt(1) = y;
}

// delete and clear the grid
void Grid::grid_clear()
{
    // clear the current cell list
    cell_list.clear();
}

// function which maps a matrix index to a cell
Cell* Grid::get_cell(const int& i, const int& j)
{
    // define helpers for number of cells
    // in each direction
    int nx = number_of_cells_x;
    int ny = number_of_cells_y;

    //return cell_list[mod(j, nx) + mod(i * nx, ny)];
    //return cell_list[mod(i, nx) + mod(j * nx, ny)];
    return cell_list[mod(i, nx) + mod(j, ny) * nx];
}

// returns the number of grid coordinates
unsigned Grid::get_ncoord()
{
    return L.size();
}

// returns the position in the grid of a particles based on
// the particles position
std::vector<int> Grid::get_cell_coordinate(Particle* particle_pt)
{
    // dereference the particle positon
    double y = *particle_pt->q_pt(1);
    double x = *particle_pt->q_pt(0) - get_translational_shift(y);

    // make a vector of coordinates
    std::vector<int> cell_coordinate(2, 0);

    // get the position of the particle in the grid
    int pos_x = (int)floor(x / L[0] * number_of_cells_x);
    int pos_y = (int)floor(y / L[2] * number_of_cells_y);

    // assign the cell coordinate
    //cell_coordinate[0] = mod(pos_x, number_of_cells_x);
    //cell_coordinate[1] = mod(pos_y, number_of_cells_y);
    cell_coordinate[0] = pos_x;
    cell_coordinate[1] = pos_y;

    // return the vector with the coordinates
    return cell_coordinate;
}

// returns the distance square between two particles
std::vector<double> Grid::get_distance_square(Particle* current_particle,
                                              Particle* neighbour_particle)
{
    // make vector for all the values
    std::vector<double> distances(3, 0.0);

    // calculate x separation
    double r_x = *neighbour_particle->q_pt(0) - *current_particle->q_pt(0);

    // enforce minimum image convention
    if(r_x > 0.5 * L[0])
        r_x -= L[0];
    else if(r_x <= -0.5 * L[0])
        r_x += L[0];

    // calculate y separation
    double r_y = *neighbour_particle->q_pt(1) - *current_particle->q_pt(1);

    // enforce minimum image convention
    if(r_y > 0.5 * L[2])
        r_y -= L[2];
    else if(r_y <= -0.5 * L[2])
        r_y += L[2];

    // set the values
    distances[0] = r_x * r_x + r_y * r_y;
    distances[1] = r_x;
    distances[2] = r_y;

    // return the square of the distance
    return distances;
}

// calculates the translational shift of the box
// as a function of y which is the position of a
// particle
double Grid::get_translational_shift(const double& y)
{
    return y * L[1] / L[2];
}

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
    // calculate the correct number of cells
    int new_number_of_cells_x = (int)floor(L[0] / (2.0 * sqrt(cut_off_sq)));
    int new_number_of_cells_y = (int)floor(L[2] / (2.0 * sqrt(cut_off_sq)));

    // boolean if the grid needs to change
    bool need_to_change_grid = false;

    // perform check in x direction
    if(new_number_of_cells_x != number_of_cells_x)
    {
        need_to_change_grid = true;
        number_of_cells_x = new_number_of_cells_x;
    }

    // perform the check in the y direction
    if(new_number_of_cells_y != number_of_cells_y)
    {
        need_to_change_grid = true;
        number_of_cells_y = new_number_of_cells_y;
    }

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
      // set the neightbours
      get_cell(i, j)->set_neighbour(0, get_cell(i+1, j-1));
      get_cell(i, j)->set_neighbour(1, get_cell(i, j));
      get_cell(i, j)->set_neighbour(2, get_cell(i+1, j));
      get_cell(i, j)->set_neighbour(3, get_cell(i, j+1));
      get_cell(i, j)->set_neighbour(4, get_cell(i+1, j+1));
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

  // set the holders for the position of the particles
  double* x = NULL;
  double* y = NULL;

  double* p_x = NULL;
  double* p_y = NULL;

  double mx = 0.0;
  double my = 0.0;

  Particle* particle_k = NULL;

  // particles per direction
  int particles_in_dir = (int)ceil(sqrt(number_of_particles));

  // grid separation
  // double separation_in_x = L[0] / (particles_in_dir + 1);
  // double separation_in_y = L[2] / (particles_in_dir + 1);

  // hexagonal separation
  double separation_in_x = 2.0 * L[0] / (3.0 *(particles_in_dir + 1.0));
  double separation_in_y = L[2] / (particles_in_dir + 1);


  // particle counter
  unsigned k = 0;

  // loop over each direction
  for(int j=0; j<particles_in_dir; j++)
    for(int i=0; i<particles_in_dir; i++)
    {
        if(k < number_of_particles)
        {
            // derefernce the particle
            particle_k = molecule_pt->particle_pt(k);

            // get pointers to the partilce position
            x = particle_k->q_pt(0);
            y = particle_k->q_pt(1);

            // set the particle position grid
            // *y = (j+1) * separation_in_y;
            // *x = (i+1) * separation_in_x + get_translational_shift(*y);

            // set particle position hexagonal
            *y = (j+1) * separation_in_y;
            if(j % 2 == 0)
                *x = ((i+1.5) + (i % 2 != 0 && i != 0) * (i)/2 + (i % 2 == 0 && i > 1) * (i+1) / 2 ) * separation_in_x + get_translational_shift(*y);
            else if (j % 2 == 1)
                *x = ((i+1) + (i % 2 != 0 && i != 0) * (i+1)/2 + (i % 2 == 0 && i > 1) * (i) / 2 ) * separation_in_x + get_translational_shift(*y);
            // else if (j % 3 == 2)
            //     *x = ((i+1) + (i % 2 != 0 && i != 0) * (i+1)/2 + (i % 2 == 0 && i > 1) * (i+1) / 2 ) * separation_in_x + get_translational_shift(*y);

            // printf("%.0d mapping to %.0d\n",i, (i+1) + (i % 2 != 0 && i != 0) * (i+1)/2 + (i % 2 == 0 && i > 1) * (i) / 2);

            // get pointers to the particle momentum
            p_x = particle_k->p_pt(0);
            p_y = particle_k->p_pt(1);

            // get the mass of the particles
            mx = *particle_k->m_pt(0);
            mx = *particle_k->m_pt(1);

            // set the momentum of the particle
            // *p_x = sqrt(mx / beta) * gsl_ran_gaussian(generator, 1.0);
            // *p_y = sqrt(my / beta) * gsl_ran_gaussian(generator, 1.0);
            *p_x = sqrt(mx / beta) * normal_momentum();
            *p_y = sqrt(my / beta) * normal_momentum();

        // if particle rotates set all the rotations
        if(particle_k->rigid_body() != false)
        {
            double* Q_00 = particle_k->Q_pt(0, 0);
            double* Q_01 = particle_k->Q_pt(0, 1);
            double* Q_10 = particle_k->Q_pt(1, 0);
            double* Q_11 = particle_k->Q_pt(1, 1);

            // momentum
            //double* pi = particle_k->pi_pt(0);

            // angle
            //double alpha = M_PI * gsl_ran_gaussian(generator, 1.0);
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
            // else if(j % 3 == 2)
            // {
            //     if(i % 2 == 0)
            //         alpha = M_PI/3.0;
            // }


            // set the rotation matrix
            *Q_00 = cos(alpha);
            *Q_01 = -sin(alpha);
            *Q_10 = sin(alpha);
            *Q_11 = cos(alpha);

            // double I = *particle_k->I_pt(0);

            // set the momentum
            //*pi = sqrt(I / beta) * gsl_ran_gaussian(generator, 1.0);
        }
      }
      // increment the particle counter
      k++;
    }
}

// update the position of all the particles on the grid
void Grid::update_particles_on_grid()
{
  // make parameters
  ListNode* cell_conductor = NULL;
  ListNode* cell_next = NULL;

  Particle* particle_k = NULL;
  std::vector<int> cell_coordinate(2,0);

  // loop over all the cells
  for(int j=0; j<number_of_cells_y; j++)
    for(int i=0; i<number_of_cells_x; i++)
    {
      // get the conductor for this cell
      cell_conductor = get_cell(i, j)->get_particle_list_head();

      // loop over all the particles in this cell
      while(cell_conductor != NULL)
      {
        // derefernce the next particle
        // (need to hold this in case we move particle such that the list is
        // not changed)
        cell_next = cell_conductor->next;

        // dereference the partilce
        particle_k = cell_conductor->particle;

        // if we are in boundary cell then enforce the
        // periodic boundary condition
        if(i == 0 or i == number_of_cells_y-1)
          enforce_periodic_particle_boundary_condition(particle_k);
        else if(j == 0 or j == number_of_cells_x-1)
          enforce_periodic_particle_boundary_condition(particle_k);

        // get the potentially new cell cordinates
        cell_coordinate = get_cell_coordinate(particle_k);

        // check if particle has moved cell
        if(cell_coordinate[0] != i or cell_coordinate[1] != j)
        {
          // dereference for readability
          int new_i = cell_coordinate[0];
          int new_j = cell_coordinate[1];

          // move the particle to the new cell
          get_cell(i,j)->move_list_node(cell_conductor,
                                        get_cell(new_i, new_j));
        }

        // step the conductor forward
        cell_conductor = cell_next;
      }
    }
}

// access for grid dimensions
double* Grid::L_pt(const unsigned& i)
{
  return &L[i];
}

// access for grid dimension momentum
double* Grid::Lp_pt(const unsigned& i)
{
  return &Lp[i];
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

    // mass and momentum for particle k
    std::vector<double> momentum(2, 0.0); // senare might break for higher dim
    double massx = 0.0;
    double massy = 0.0;

    double pi = 0.0;
    double I = 0.0;

    // helpers
    Particle* particle_k = NULL;
    ListNode* cell_conductor = NULL;

    double N = 1.0;

    // loop over all the cells
    for(int j=0; j<number_of_cells_y; j++)
      for(int i=0; i<number_of_cells_x; i++)
      {
          // get the conductor for this cell
          cell_conductor = get_cell(i, j)->get_particle_list_head();

          // loop over all the particles in this cell
          while(cell_conductor != NULL)
          {
              // dereference the partilce
              particle_k = cell_conductor->particle;

              //momentum = get_box_momentum(particle_k);
              momentum[0] = *particle_k->p_pt(0);
              momentum[1] = *particle_k->p_pt(1);

              // mass in both directions#pragma omp for
              massx = *particle_k->m_pt(0);
              massy = *particle_k->m_pt(1);

              // add to the temperature
              temperature = momentum[0] * momentum[0] / (massx * N)
                          + momentum[1] * momentum[1] / (massy * N)
                          + (N-1.0)/N * temperature;

              //if particle rotates add to rotational temperature
              if(particle_k->rigid_body())
              {
                  pi = *particle_k->pi_pt(0);
                  I =  *particle_k->I_pt(0);

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
        total_temp = 0.5 * temperature;

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
  // vector for holding distance
  std::vector<double> distance(3, 0.0);

  // calculate the square of the distance between particles
  distance = get_distance_square(current_particle, neighbour_particle);

  // chcek the cutoff criterion
  if(distance[0] < cut_off_sq)
  {
    // calculate the actual distance
    distance[0] = sqrt(distance[0]);

    // Calculate force and potential
    system_pt->compute_pair_force(molecule_pt, current_particle,
                                  neighbour_particle, distance[0], distance[1],
                                  distance[2]);
  }
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

  // parameters for cell loops
  Cell* current_cell = NULL;
  Cell* neighbour_cell = NULL;

  ListNode* current_conductor = NULL;
  ListNode* neighbour_conductor = NULL;

  Particle* current_particle = NULL;
  Particle* neighbour_particle = NULL;

  // make iterators to do Newton iteration
  // forces are equal and opposite
  unsigned current_newton_iterator = 0;
  unsigned neighbour_newton_iterator = 0;

  // loop over all the boxes to calculate the forces
#pragma omp for
   for(int j=0; j<number_of_cells_y; j++)
   {
     for(int i=0; i<number_of_cells_x; i++)
     {
       // get the current cell
       current_cell = get_cell(i,j);

       for(unsigned n=0; n<number_of_neighbours; n++)
       {
         // get the current cells conductor
         current_conductor = current_cell->get_particle_list_head();

         // get the neighbour cell
         neighbour_cell = current_cell->get_neighbour(n);

         // reset newton iterator
         current_newton_iterator = 0;

         // loop over particles in current cell
         while(current_conductor != NULL)
         {
           // conductor over neighbour cell
           neighbour_conductor = neighbour_cell->get_particle_list_head();

           // get the particle pointer
           current_particle = current_conductor->particle;

           // reset newton iterator
           neighbour_newton_iterator = 0;

           // loop over particles in neighbour cell
           while(neighbour_conductor != NULL)
           {
             // apply Newton iterators if we are i the same cell
             if(current_cell == neighbour_cell)
             {
               if(current_newton_iterator <= neighbour_newton_iterator)
               {
                  // get the neighbour particle
                  neighbour_particle = neighbour_conductor->particle;

                  // don't compute the force for the same particles
                  if(current_particle != neighbour_particle)
                    compute_force(system_pt, molecule_pt,
                                  current_particle, neighbour_particle);
                }

                // iterate the newton current newton iterator
                neighbour_newton_iterator++;
             }
             else
             {
               // get the neighbour particle
               neighbour_particle = neighbour_conductor->particle;

               // compute the force
               compute_force(system_pt, molecule_pt, current_particle,
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
               L[i] = parsed_row[i];
           else
               Lp[i-3] = parsed_row[i];
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

   Particle* particle_k = NULL;
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

      particle_k = molecule_pt->particle_pt(k);

      // set position
      *particle_k->q_pt(0) = parsed_row[0];
      *particle_k->q_pt(1) = parsed_row[1];

      // if we rotates
      if(particle_k->rigid_body())
      {
          *particle_k->Q_pt(0, 0) = parsed_row[2];
          *particle_k->Q_pt(0, 1) = parsed_row[3];

          *particle_k->Q_pt(1, 0) = parsed_row[4];
          *particle_k->Q_pt(1, 1) = parsed_row[5];
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

    Particle* particle_k = NULL;
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

       particle_k = molecule_pt->particle_pt(k);

       // set momentum
       *particle_k->p_pt(0) = parsed_row[0];
       *particle_k->p_pt(1) = parsed_row[1];

       // if we rotates
       if(particle_k->rigid_body())
       {
           *particle_k->pi_pt(0) = parsed_row[2];
       }

       // increment particle counter
       k++;
      }
}
