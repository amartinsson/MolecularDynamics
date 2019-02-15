#include "NptGrid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                               NPT GRID CLASS
// ------------------------------------------------------------------------- //
NptGrid::NptGrid(const double& a_x, const double& b_x, const double& b_y,
                 const double& cut_off, Molecule* molecule_pt,
                 const double& mass, const double& target_press) :
                 Grid(a_x, b_x, b_y, cut_off, molecule_pt)
{
    // set up tracking object to track the average
    // of the temperature
    Volume_pt = new AverageObservable();
    Pressure_pt = new AverageObservable();

    // set the mass of the cell grid for the
    // NPT integration
    box_mass = mass;

    // set the target pressure of the cell grid
    // for comupting the box force
    target_pressure = target_press;

    // reserve space for the box force in each direction
    box_grad_potential.resize(3);
    //box_grad_potential = new double[3];
    momentum_sq.resize(3);

    for(unsigned i=0; i<momentum_sq.size(); i++)
    {
        box_grad_potential[i] = 0.0;
        momentum_sq[i] = 0.0;
    }

    printf("Inital grid set as %d x %d\n",
           number_of_cells_x, number_of_cells_y);
}

// destructor
NptGrid::~NptGrid()
{
    // delete all the average observables
    delete Volume_pt;
    delete Pressure_pt;

    delete &box_grad_zero;
    delete &box_grad_one;
    delete &box_grad_two;

    // remove all the box force variables
    for(unsigned i=0; i<momentum_sq.size(); i++)
    {
        delete &box_grad_potential[i];
        delete &momentum_sq[i];
    }

    // clear all the vectors
    box_grad_potential.clear();
    momentum_sq.clear();
}

// computes the force in and respects the cut off radius. This is
// a helper function because the iteration in the same cell is different
// as we should not double count some distances
void NptGrid::compute_force(System* system_pt, Molecule* molecule_pt,
                            Particle* current_particle,
                            Particle* neighbour_particle)
{
    // int test = 0;
    // #pragma omp target data map(to:test)
    // {
    //     printf("This is GPU\n");
    //     while(test != 2)
    //     {
    //         if(test == 0)
    //             test += 1;
    //         else
    //             test -= 1;
    //     }
    // }
  // vector for holding distance
  std::vector<double> distance(3, 0.0);

  // calculate the square of the distance between particles
  distance = get_distance_square(current_particle, neighbour_particle);

  // chcek the cutoff criterion
  if(distance[0] < cut_off_sq)
  {
      // holder for the pair force
      vector<double> f_ij(2, 0.0);
      vector<double> r_tilde(2, 0.0);

      // calculate the actual distance
      distance[0] = sqrt(distance[0]);

      // Calculate force and potential
      f_ij = system_pt->compute_pair_force(molecule_pt, current_particle,
                                           neighbour_particle, distance[0],
                                           distance[1], distance[2]);

      // rescale the separation into box invariant coordinates
      r_tilde = get_box_min_image_sep(current_particle, neighbour_particle);

      // these need to be minus additive as we are calculating the
      // force but we need the gradient!
      //
      // These are the virials
      box_grad_zero += r_tilde[0] * f_ij[0];
      box_grad_one  += r_tilde[1] * f_ij[0];
      box_grad_two  += r_tilde[1] * f_ij[1];
  }
}

  // function which calculates and returns the
  // volume of the current cell
  double NptGrid::get_instant_volume()
  {
    return Volume_pt->get_instant();
  }

  // function which calculates and returns the
  // volume of the current cell
  double NptGrid::get_instant_pressure()
  {
    return Pressure_pt->get_instant();
  }

  // return the volume
  double NptGrid::get_volume()
  {
    return Volume_pt->get_average();
  };

  // returns the currently held pressure
  double NptGrid::get_pressure()
  {
    return Pressure_pt->get_average();
  }

  // returns the mass of the cell grid
  double NptGrid::get_mass()
  {
    return box_mass;
  }

  // function which updates the pressure
  void NptGrid::update_pressure()
  {
    // calculate the macroscopic pressure
    double pressure = calculate_pressure();

    // add to the observables
    Pressure_pt->observe(pressure);
  }

  // enforces the relaative positoon of the particles
  void NptGrid::enforce_constant_relative_particle_pos(
                                            const std::vector<double>& L_old)
  {
      // make sure particle is in box
      update_particles_on_grid();

      // rescale positions
      enforce_relative_particle(L_old);
  }

  std::vector<double> NptGrid::get_box_min_image_sep(Particle* current_particle,
                                                   Particle* neighbour_particle)
  {
      // get the box coordinates
      std::vector<double> qc = get_box_coordinate(current_particle);
      std::vector<double> qn = get_box_coordinate(neighbour_particle);

      // calculate x separation
      double r_x = qn[0] - qc[0];

      // enforce minimum image convention
      if(r_x > 0.5)
        r_x -= 1.0;
      else if(r_x <= -0.5)
        r_x += 1.0;

      // calculate y separation
      double r_y = qn[1] - qc[1];

      // enforce minimum image convention
      if(r_y > 0.5)
        r_y -= 1.0;
      else if(r_y <= -0.5)
        r_y += 1.0;

      // make new vector
      std::vector<double> r_tilde(2, 0.0);

      // assign the differeces
      r_tilde[0] = r_x;
      r_tilde[1] = r_y;

      // retrun the rescalied vector
      return r_tilde;
  }

  // This function returns the non dimensional coordinate w.r.t the box
  // of a particle in the box.
  std::vector<double> NptGrid::get_box_coordinate(Particle* particle_pt)
  {
    // make coordinate vector
    std::vector<double> q_tilde(2, 0.0);

    // get the particle coordinate
    double x = *particle_pt->q_pt(0);
    double y = *particle_pt->q_pt(1);

    // get box coordiantes
    double a_x = *L_pt(0);
    double b_x = *L_pt(1);
    double b_y = *L_pt(2);

    // calculate the coordinates
    q_tilde[0] = (1.0 / a_x) * x - (b_x / (a_x * b_y)) * y;
    q_tilde[1] = (1.0 / b_y) * y;

    // return the coordinates
    return q_tilde;
  }

  // function which sets the position of a particle to the invariant
  // position given by q tilde
  void NptGrid::set_box_coordinate(std::vector<double>& q_tilde,
                                   Particle* particle_pt)
  {
      // get the box size
      double ax = *L_pt(0);
      double bx = *L_pt(1);
      double by = *L_pt(2);

      double* x = particle_pt->q_pt(0);
      double* y = particle_pt->q_pt(1);

      // set both directions
      *x = ax * q_tilde[0] + bx * q_tilde[1];
      *y = by * q_tilde[1];
  }

  // This function calculates the non dimensional momentum coordinate w.r.t
  // to the box of the particle given
  std::vector<double> NptGrid::get_box_momentum(Particle* particle_pt)
  {
    // make coordinate vector
    std::vector<double> p_tilde(2, 0.0);

    // get the particle coordinate
    double p_x = *particle_pt->p_pt(0);
    double p_y = *particle_pt->p_pt(1);

    // get box coordiantes
    double a_x = *L_pt(0);
    double b_x = *L_pt(1);
    double b_y = *L_pt(2);

    // calculate the scaled momentum
    p_tilde[0] = a_x * p_x + b_x * p_y;
    p_tilde[1] = b_y * p_y;

    // return the scaled momentum
    return p_tilde;
  }

  // function which sets the position of a particle to the invariant
  // position given by q tilde
  void NptGrid::set_box_momentum(std::vector<double>& p_tilde,
                                 Particle* particle_pt)
  {
      // get the box size
      double ax = *L_pt(0);
      double bx = *L_pt(1);
      double by = *L_pt(2);

      double* px = particle_pt->p_pt(0);
      double* py = particle_pt->p_pt(1);

      // set both directions
      *px = p_tilde[0]/ax - bx/(ax*by) * p_tilde[1];
      *py = p_tilde[1]/by;
  }

  // function which update the instantenous and average volume
  void NptGrid::update_volume()
  {
    // derefernce the cell
    double a_x = *L_pt(0);
    double b_y = *L_pt(2);

    // calcuate the volume
    double volume = a_x * b_y;

    // make observation of the volume
    Volume_pt->observe(volume);
  }

  // update both the pressure and the temperature variables -
  // this needs to be called at the end of every integration step
  // and will update the pressure and temperature
  void NptGrid::update_pressure_temperature()
  {
    // first update the volume
    update_volume();

    // second update the temperature
    update_temperature();

    // update the pressure assuming that the q_dot_f value has been updated
    update_pressure();
  }


  // function which calculates the pressure
  double NptGrid::calculate_pressure()
  {
      // dereference the box
      double lx1 = *L_pt(0);
      double ly1 = *L_pt(1);
      double ly2 = *L_pt(2);

      double pressure = 1.0/(lx1*ly2) * momentum_sq[0]
                        - 1.0/(ly2) * box_grad_potential[0]
                        + 1.0/(lx1*ly2) * momentum_sq[2]
                        - ly1/(ly2*lx1*lx1) * momentum_sq[1]
                        - 1.0/(lx1) * box_grad_potential[2];

      // return the calulcated pressure
      return 0.5 * pressure;
  }

  // enforce the constraint that the relative distances
  // of the particle position and momentum cannot change.
  void NptGrid::enforce_relative_particle(const std::vector<double> L_old)
  {
      // make parameters
      ListNode* cell_conductor = NULL;
      Particle* particle_k = NULL;

      // get box coordiantes
      double ax = *L_pt(0);
      double bx = *L_pt(1);
      double by = *L_pt(2);

      // dereference the old box coordinate
      double ax_old = L_old[0];
      double bx_old = L_old[1];
      double by_old = L_old[2];

      // reference for particles data
      double* qx = NULL;
      double* qy = NULL;

      double* px = NULL;
      double* py = NULL;

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

            // dereference the position and momentum
            qx = particle_k->q_pt(0);
            qy = particle_k->q_pt(1);

            px = particle_k->p_pt(0);
            py = particle_k->p_pt(1);

            // set positions
            *qx = ax/ax_old * (*qx)
                  + (bx/by_old - ax * bx_old/(ax_old * by_old)) * (*qy);

            *qy = by/by_old * (*qy);

            // set momentums
            *px = ax_old/ax * (*px)
                  + (bx_old/ax - bx * by_old/(ax * by)) * (*py);

            *py = by_old/by * (*py);

            // step the conductor forward
            cell_conductor = cell_conductor->next;
          }
        }
  }

  // updates the box force in all directions
  void NptGrid::update_accumulted_momentum()
  {
    // reset the values of the vectors
    for(unsigned i=0; i< momentum_sq.size(); i++)
    {
        momentum_sq[i] = 0.0;
    }

    // reference pointers
    ListNode* cell_conductor = NULL;
    Particle* particle_k = NULL;

    // helpers
    double mx = 0.0;
    double my = 0.0;

    double px = 0.0;
    double py = 0.0;

    // loop over all the particles
    for(int j=0; j<number_of_cells_y; j++)
      for(int i=0; i<number_of_cells_x; i++)
      {
        // get the conductor for this cell
        cell_conductor = get_cell(i, j)->get_particle_list_head();

        // loop over all particles in the cell
        while(cell_conductor != NULL)
        {
          // dereference the partilce
          particle_k = cell_conductor->particle;

          // dereference the particle
          mx = *particle_k->m_pt(0);
          my = *particle_k->m_pt(1);

          px = *particle_k->p_pt(0);
          py = *particle_k->p_pt(1);

          // update the values
          momentum_sq[0] += px * px / mx;
          momentum_sq[1] += px * py / mx;
          momentum_sq[2] += py * py / my;

          // increment the conductor
          cell_conductor = cell_conductor->next;
        }
      }
  }

  // fuction which updates the particle forces. It automatically checks if the
  // grid has changed sufficently to have to rebuild the grid. It then contineous
  // updateing the forces of all the particles using the linked lists from above,
  // looping over all the cells in the grid. It also updates the pressure and
  // the temperature by calculating the dot product of q and f for all the
  // particles
  void NptGrid::update_particle_forces(System* system_pt,
                                       Molecule* molecule_pt)
  {
    // if the grid needs to be rebuilt do that.
    // Otherwise simply make sure that all the particles
    // are in the correct cells and that the boundary
    // conditions are enforced.
  // #pragma omp single
    {
      if(rebuild_grid_check())
      {
        rebuild_periodic_grid(molecule_pt);
        printf("\tfinished rebuilding %d x %d grid\n",
               number_of_cells_x, number_of_cells_y);
      }
      else
        update_particles_on_grid();

      // clear all the particle forces
      clear_particle_forces(molecule_pt);

      // reset the potential gradient w.r.t to the box coordinate
      box_grad_potential[0] = 0.0;
      box_grad_potential[1] = 0.0;
      box_grad_potential[2] = 0.0;
    }

    box_grad_zero = 0.0;
    box_grad_one = 0.0;
    box_grad_two = 0.0;

    // loop over all the boxes to calculate the forces
#pragma omp parallel for \
        reduction(+:box_grad_zero, box_grad_one, box_grad_two) \
        schedule(dynamic) collapse(3) \
        firstprivate(number_of_cells_x, number_of_neighbours, system_pt)
     for(int j=0; j<number_of_cells_y; j++)
     {
       for(int i=0; i<number_of_cells_x; i++)
       {
         for(unsigned n=0; n<number_of_neighbours; n++)
         {
             // get the current cell
             Cell* current_cell = get_cell(i,j);

             // get the current cells conductor
             ListNode* current_conductor
                                    = current_cell->get_particle_list_head();

             // get the neighbour cell
             Cell* neighbour_cell = current_cell->get_neighbour(n);

             // reset newton iterator
             unsigned current_newton_iterator = 0;

             // loop over particles in current cell
             while(current_conductor != NULL)
             {
                 // conductor over neighbour cell
                 ListNode* neighbour_conductor
                                    = neighbour_cell->get_particle_list_head();

                 // get the particle pointer
                 Particle* current_particle = current_conductor->particle;

                 // reset newton iterator
                 unsigned neighbour_newton_iterator = 0;

                 // loop over particles in neighbour cell
                 while(neighbour_conductor != NULL)
                 {
                     // apply Newton iterators if we are in the same cell
                     if(current_cell == neighbour_cell)
                     {
                         if(current_newton_iterator <= neighbour_newton_iterator)
                         {
                             // get the neighbour particle
                             Particle* neighbour_particle
                                            = neighbour_conductor->particle;

                             // don't compute the force for the same particles
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
                         Particle* neighbour_particle
                                            = neighbour_conductor->particle;

                         // compute the force
                         compute_force(system_pt, molecule_pt,
                                       current_particle, neighbour_particle);
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

 // set the reduction variables
 box_grad_potential[0] = box_grad_zero;
 box_grad_potential[1] = box_grad_one;
 box_grad_potential[2] = box_grad_two;
}
