#include "NptGrid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                               NPT GRID CLASS
// ------------------------------------------------------------------------- //
NptGrid::NptGrid(const double& a_x, const double& b_x, const double& b_y,
                 const double& cut_off, Molecule* molecule_pt,
                 const double& mass, const double& target_press) :
                 Grid(a_x, b_x, b_y, cut_off, molecule_pt), momentum_sq(3),
                 virial(2,2), nablaSa(2,2), nablaSb(2,2), nablaSd(2,2)
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
    // box_grad_potential.resize(3);
    //box_grad_potential = new double[3];
    // momentum_sq.resize(3);

    // for(unsigned i=0; i<momentum_sq.size(); i++)
    //     box_grad_potential[i] = 0.0;

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
  distance = get_distance_square(*current_particle, *neighbour_particle);

  // chcek the cutoff criterion
  if(distance[0] < cut_off_sq)
  {
      // holder for the pair force
      vector<double> f_ij(2, 0.0);
      Vector r_tilde(2);

      // calculate the actual distance
      distance[0] = sqrt(distance[0]);

      // Calculate force and potential
      f_ij = system_pt->compute_pair_force(molecule_pt, current_particle,
                                           neighbour_particle, distance[0],
                                           distance[1], distance[2]);

      // rescale the separation into box invariant coordinates
      r_tilde = get_box_min_image_sep(*current_particle, *neighbour_particle);

      // these need to be minus additive as we are calculating the
      // force but we need the gradient!
      //
      // These are the virials
      // box_grad_zero += r_tilde(0) * f_ij[0];
      // box_grad_one  += r_tilde(1) * f_ij[0];
      // box_grad_two  += r_tilde(1) * f_ij[1];
      box_grad_zero -= r_tilde(0) * f_ij[0];
      box_grad_one  -= r_tilde(1) * f_ij[0];
      box_grad_two  -= r_tilde(1) * f_ij[1];
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
      // updatae the accumulated momentum
      update_accumulted_momentum();

      // calculate the macroscopic pressure
      double pressure = calculate_pressure();

      // add to the observables
      Pressure_pt->observe(pressure);
  }

  // enforces the relaative positoon of the particles
  void NptGrid::enforce_constant_relative_particle_pos(const Matrix& Sold)
  {
      // make sure particle is in box
      update_particles_on_grid();

      // rescale positions
      enforce_relative_particle(Sold);
  }

Vector NptGrid::get_box_min_image_sep(const Particle& current_particle,
                                      const Particle& neighbour_particle)
  {
      // get the box coordinates
      Vector qc = get_box_coordinate(current_particle);
      Vector qn = get_box_coordinate(neighbour_particle);

      // calculate x separation
      Vector r = qn - qc;

      // enforce minimum image convention
      calculate_min_image(r);

      // retrun the rescalied vector
      return r;
  }

  // This function returns the non dimensional coordinate w.r.t the box
  // of a particle in the box.
  Vector NptGrid::get_box_coordinate(const Particle& particle)
  {
    // make coordinate vector
    // Vector q_tilde(2);

    // rescale the coordiinate
    Vector q_tilde = S.inv() * particle.q;

    // return the coordinates
    return q_tilde;
  }

  // function which sets the position of a particle to the invariant
  // position given by q tilde
  void NptGrid::set_box_coordinate(const Vector& q_tilde, Particle& particle)
  {
      // calculate new position
      particle.q = S * q_tilde;
  }

  // This function calculates the non dimensional momentum coordinate w.r.t
  // to the box of the particle give
Vector NptGrid::get_box_momentum(const Particle& particle)
{
    // calulate the new momentum
    Vector p_tilde = S * particle.p;

    // return the scaled momentum
    return p_tilde;
}

  // function which sets the position of a particle to the invariant
  // position given by q tilde
  void NptGrid::set_box_momentum(const Vector& p_tilde, Particle& particle)
  {
      // set the new momentum
      particle.p = S.inv() * p_tilde;
      //
      // // get the box size
      // double ax = *L_pt(0);
      // double bx = *L_pt(1);
      // double by = *L_pt(2);
      //
      // double* px = particle_pt->p_pt(0);
      // double* py = particle_pt->p_pt(1);
      //
      // // set both directions
      // *px = p_tilde[0]/ax - bx/(ax*by) * p_tilde[1];
      // *py = p_tilde[1]/by;
  }

  // function which update the instantenous and average volume
  void NptGrid::update_volume()
  {
    // derefernce the cell
    // double a_x = *L_pt(0);
    // double b_y = *L_pt(2);

    // calcuate the volume
    // double volume = a_x * b_y;
    double volume = S.det();

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

  // this function update the gradient matrices that are used for finding the
  // momentum of the box
  void NptGrid::update_gradient_matrices()
  {
      // update all three matrices
      nablaSa(0,0) = -1.0 / Grid::S(0,0);

      nablaSb(0,1) = - 0.5 / Grid::S(0,0);
      nablaSb(1,0) = - 0.5 / Grid::S(0,0);

      nablaSd(0,1) = 0.5 * Grid::S(0,1) / (Grid::S(0,0) * Grid::S(1,1));
      nablaSd(1,0) = 0.5 * Grid::S(0,1) / (Grid::S(0,0) * Grid::S(1,1));
      nablaSd(1,1) = - 1.0 / Grid::S(1,1);
  }

  // function which calculates the pressure
  double NptGrid::calculate_pressure()
  {
      // // dereference the box
      // double lx1 = S(0, 0);
      // double ly1 = S(0, 1);
      // double ly2 = S(1, 1);
      //
      // double pressure = 1.0/(lx1*ly2) * momentum_sq(0)
      //                   - 1.0/(ly2) * virial(0,0)
      //                   + 1.0/(lx1*ly2) * momentum_sq(2)
      //                   - ly1/(ly2*lx1*lx1) * momentum_sq(1)
      //                   - 1.0/(lx1) * virial(1,1);
      //
      // // return the calulcated pressure
      // return 0.5 * pressure;

      // double pressure = -1.0/S(1,1) * (momentum_sq(0) + virial(0,0) + target_pressure * S(1,1))
      //                    -1.0/S(0,0) * (momentum_sq(1) + virial(1,1) + target_pressure * S(0,0));
      // return pressure;
    double pressure = -0.5/S(1,1) * (momentum_sq(0) + virial(0,0))
                      -0.5/S(0,0) * (momentum_sq(1) + virial(1,1));
    if(pressure > 10.0 * target_pressure)
    {
        printf("\n ERROR: Target pressure is: %1.3f breaking!\n\n", pressure);
        exit(-1);
    }

    return pressure;
  }

// enforce the constraint that the relative distances
// of the particle position and momentum cannot change.
void NptGrid::enforce_relative_particle(const Matrix& Sold)
{
    Matrix qScale(2,2);
    Matrix pScale(2,2);

    qScale = S * Sold.inv();
    pScale = S.inv() * Sold;

    // make parameters
    ListNode* cell_conductor = NULL;
    Particle* particle_k = NULL;

    // make holders
    Vector qnew(2);
    Vector pnew(2);

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

                 qnew = qScale * (*particle_k).q;
                 pnew = pScale * (*particle_k).p;

                 (*particle_k).q = qnew;
                 (*particle_k).p = pnew;

                 // step the conductor forward
                 cell_conductor = cell_conductor->next;
            }
        }
  }

// updates the box force in all directions
void NptGrid::update_accumulted_momentum()
{
    // reset the values of the vectors
    for(unsigned i=0; i<momentum_sq.size(); i++)
        momentum_sq(i) = 0.0;

    // reference pointers
    ListNode* cell_conductor = NULL;
    Particle* particle;

    // loop over all the particles
#pragma omp simd collapse(2)
    for(int j=0; j<number_of_cells_y; j++)
        for(int i=0; i<number_of_cells_x; i++)
        {
            // get the conductor for this cell
            cell_conductor = get_cell(i, j)->get_particle_list_head();

            // loop over all particles in the cell
            while(cell_conductor != NULL)
            {
                // dereference the partilce
                particle = cell_conductor->particle;

                Matrix Ka = this->nablaSa / particle->m(0,0);
                Matrix Kd = this->nablaSd;

                // divide by the correct mass
                Kd(0,1) /= particle->m(0,0);
                Kd(1,0) /= particle->m(0,0);
                Kd(1,1) /= particle->m(1,1);

                momentum_sq(0) += (particle->p.T() * Ka).dot(particle->p);
                momentum_sq(1) += (particle->p.T() * Kd).dot(particle->p);

                // momentum_sq(0) += particle->p(0) * particle->p(0)
                //                 / particle->m(0,0);
                // momentum_sq(1) += particle->p(0) * particle->p(1)
                //                 / particle->m(0,0);
                // momentum_sq(2) += particle->p(1) * particle->p(1)
                //                 / particle->m(1,1);

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
        // printf("\tfinished rebuilding %d x %d grid\n",
        //        number_of_cells_x, number_of_cells_y);
      }
      else
        update_particles_on_grid();

      // clear all the particle forces
      clear_particle_forces(molecule_pt);

      // reset the potential gradient w.r.t to the box coordinate
      virial.zero();
    }

    box_grad_zero = 0.0;
    box_grad_one = 0.0;
    box_grad_two = 0.0;

    // loop over all the boxes to calculate the forces
#pragma omp parallel for default(shared)\
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
 virial(0,0) = box_grad_zero;
 virial(0,1) = box_grad_one;
 virial(1,1) = box_grad_two;
}
