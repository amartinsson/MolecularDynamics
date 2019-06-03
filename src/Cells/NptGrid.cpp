#include "NptGrid.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                               NPT GRID CLASS
// ------------------------------------------------------------------------- //
NptGrid::NptGrid(const Matrix& Szero, const double& cut_off,
                 Molecule* molecule_pt, const double& mass,
                 const double& target_press) :
                 Grid(Szero, cut_off, molecule_pt), nablaK(2,2), virial(2,2)
{
    if(Szero.size()[0] > 2)
    {
        virial.resize(3, 3);
        nablaK.resize(3, 3);
    }

    // set up tracking object to track the average of the temperature
    Volume_pt = new AverageObservable();
    Pressure_pt = new AverageObservable();

    // set the mass of the cell grid for the NPT integration
    box_mass = mass;

    // set the target pressure of the cell grid for comupting the box force
    target_pressure = target_press;
}

// destructor
NptGrid::~NptGrid()
{
    // delete all the average observables
    delete Volume_pt;
    delete Pressure_pt;
}

// computes the force in and respects the cut off radius. This is
// a helper function because the iteration in the same cell is different
// as we should not double count some distances
void NptGrid::compute_force(System* system_pt, Molecule* molecule_pt,
                            Particle* current_particle,
                            Particle* neighbour_particle)
{
    // calculate the square of the distance between particles
    Vector r = get_separation(*current_particle, *neighbour_particle);

    double rsq = r.l22();

    // chcek the cutoff criterion
    if(rsq < cut_off_sq)
    {
        // Calculate force and potential
        Vector f_ij = system_pt->compute_force(molecule_pt,
                                               current_particle,
                                               neighbour_particle,
                                               std::sqrt(rsq), r);

       //  Vector q1 = S.inv() * (*current_particle).q;
       //  Vector q2 = S.inv() * (*neighbour_particle).q;
       //
       // if(std::sqrt(rsq) < pow(2.0,1.0/6.0))
       //      printf("r = %1.3e f = %1.3e %p = (%1.3f, %1.3f, %1.3f) %p = (%1.3f, %1.3f, %1.3f)\n",
       //      std::sqrt(rsq), f_ij.l2(), current_particle, q1(0), q1(1), q1(2), neighbour_particle, q2(0), q2(1), q2(2));


        if(Grid::with_radial_dist)
            Grid::Radial_pt->update(std::sqrt(rsq));

        // rescale the separation into box invariant coordinates
        Vector r_tilde = S.inv() * r;

        // calculate the virial
        // Matrix V = f_ij.out(r_tilde);
        Matrix V(3,3);

        for(unsigned i=0; i<3; i++)
            for(unsigned j=i; j<3; j++)
            {
                Matrix K(3,3);
                K(i,j) = 1.0;
                V(i,j) = f_ij.dot(K * r_tilde);
            }

        // these need to be minus additive as we are calculating the
        // force but we need the gradient!
        //
        // These are the virials
        box_grad_00 -= V(0,0);
        box_grad_01 -= V(0,1);
        box_grad_11 -= V(1,1);

        if(r_tilde.size() > 2)
        {
            box_grad_02 -= V(0,2);
            box_grad_12 -= V(1,2);
            box_grad_22 -= V(2,2);
        }
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
  void NptGrid::enforce_constant_relative_particle_pos(const Matrix& Sold)
  {
      // make sure particle is in box
      update_particles_on_grid();

      // rescale positions
      enforce_relative_particle(Sold);
  }

// Vector NptGrid::get_box_min_image_sep(const Particle& current_particle,
//                                       const Particle& neighbour_particle)
//   {
//       // get the box coordinates
//       Vector qc = get_box_coordinate(current_particle);
//       Vector qn = get_box_coordinate(neighbour_particle);
//
//       // calculate x separation
//       Vector r = qn - qc;
//
//       // enforce minimum image convention
//       calculate_min_image(r);
//
//       // retrun the rescalied vector
//       return r;
//   }

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
    Vector p_tilde = S.T() * particle.p;

    // return the scaled momentum
    return p_tilde;
}

  // function which sets the position of a particle to the invariant
  // position given by q tilde
  void NptGrid::set_box_momentum(const Vector& p_tilde, Particle& particle)
  {
      // set the new momentum
      particle.p = S.inv().T() * p_tilde;
  }

  // function which update the instantenous and average volume
  void NptGrid::update_volume()
  {
      // calculate volume
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

// function which calculates the pressure
double NptGrid::calculate_pressure()
{
    this->update_kinetic_gradient();

    double pressure = 0.0;

    if(number_of_cells_z != 0)
        pressure = -1.0 / (3.0 * S(1,1) * S(2,2)) * (nablaK(0,0) + virial(0,0))
                   -1.0 / (3.0 * S(0,0) * S(2,2)) * (nablaK(1,1) + virial(1,1))
                   -1.0 / (3.0 * S(0,0) * S(1,1)) * (nablaK(2,2) + virial(2,2));
    else
        pressure = -1.0 / (2.0 * S(1,1)) * (nablaK(0,0) + virial(0,0))
                   -1.0 / (2.0 * S(0,0)) * (nablaK(1,1) + virial(1,1));

    return pressure;
}

// enforce the constraint that the relative distances
// of the particle position and momentum cannot change.
void NptGrid::enforce_relative_particle(const Matrix& Sold)
{
    // calculate the scale matrices
    Matrix qScale = S * Sold.inv();
    Matrix pScale = S.inv().T() * Sold.T();

    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // loop over all the cells
    //#pragma omp simd collapse(3)
    #pragma omp parallel for default(shared)\
            firstprivate(zend, number_of_cells_y, number_of_cells_x, \
                         qScale, pScale) \
            schedule(dynamic) collapse(3)
    for(int k=0; k<zend; k++)
        for(int j=0; j<number_of_cells_y; j++)
            for(int i=0; i<number_of_cells_x; i++)
                for(ListNode* cell_conductor=get_cell(i,j,k)->get_particle_list_head();
                        cell_conductor != NULL;
                            cell_conductor=cell_conductor->next)
                {
                    // dereference the partilce
                    Particle* particle = cell_conductor->particle;

                    (*particle).q = qScale * (*particle).q;
                    (*particle).p = pScale * (*particle).p;
                }
}

// update the gradient of the kinetic energy w.r.t the box
void NptGrid::update_kinetic_gradient()
{
    // clear the matrix
    nablaK.zero();
    Matrix Sinv = S.inv();

    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // loop over all the cells
    #pragma omp parallel for default(shared)\
            firstprivate(zend, number_of_cells_y, number_of_cells_x, Sinv) \
            schedule(dynamic) collapse(3)
    for(int k=0; k<zend; k++)
        for(int j=0; j<number_of_cells_y; j++)
            for(int i=0; i<number_of_cells_x; i++)
                for(ListNode* cell_conductor=get_cell(i,j,k)->get_particle_list_head();
                        cell_conductor != NULL;
                            cell_conductor=cell_conductor->next)
                {
                    // dereference the partilce
                    Particle* particle = cell_conductor->particle;

                    Matrix B = Sinv * (*particle).m.inv();

                    int dim = B.size()[0];

                    for(int m=0; m<dim; m++)
                        for(int n=m; n<dim; n++)
                        {
                            Matrix D(dim, dim);
                            D(n,m) = 1.0;
                            Matrix E = D * B;
                            #pragma omp atomic
                            nablaK(m,n) -= 0.5 * ((*particle).p.T() * (E + E.T())).dot((*particle).p);
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
    if(rebuild_grid_check())
        rebuild_periodic_grid(molecule_pt);
    else
        update_particles_on_grid();

    // clear all the particle forces
    clear_particle_forces(molecule_pt);

    // reset the potential gradient w.r.t to the box coordinate
    virial.zero();

    box_grad_00 = 0.0;
    box_grad_01 = 0.0;
    box_grad_11 = 0.0;

    box_grad_02 = 0.0;
    box_grad_12 = 0.0;
    box_grad_22 = 0.0;

    // calculate the z direction component
    int zend = 1;
    if(number_of_cells_z != 0)
        zend = number_of_cells_z;

    // loop over all the boxes to calculate the forces
#pragma omp parallel for default(shared)\
        reduction(+:box_grad_00, box_grad_01, box_grad_11, box_grad_02,\
                    box_grad_12, box_grad_22) \
        firstprivate(zend, number_of_cells_y, number_of_cells_x, system_pt) \
        schedule(dynamic) collapse(3)
    for(int k=0; k<zend; k++)
    {
        for(int j=0; j<number_of_cells_y; j++)
        {
            for(int i=0; i<number_of_cells_x; i++)
            {
                // get the current cell
                Cell* current_cell = get_cell(i, j, k);

                // for(unsigned n=0; n<number_of_neighbours; n++)
                for(const auto& ncell : current_cell->neighbour_list)
                {
                    // get the current cells conductor
                    // ListNode* current_conductor =
                                    // current_cell->get_particle_list_head();

                    // get the neighbour cell
                    Cell* neighbour_cell = ncell.second;

                    // reset newton iterator
                    unsigned current_newton_iterator = 0;

                    // loop over particles in current cell
                    // while(current_conductor != NULL)
                    // loop over particles in current cell
                    for(ListNode* current_conductor=current_cell->get_particle_list_head();
                            current_conductor != NULL;
                                current_conductor=current_conductor->next)
                    {
                        // conductor over neighbour cell
                        // ListNode* neighbour_conductor =
                        //             neighbour_cell->get_particle_list_head();

                        // get the particle pointer
                        Particle* current_particle =
                                            current_conductor->particle;

                        // reset newton iterator
                        unsigned neighbour_newton_iterator = 0;

                        // loop over particles in neighbour cell
                        // while(neighbour_conductor != NULL)
                        // loop over particles in neighbour cell
                        for(ListNode* neighbour_conductor=neighbour_cell->get_particle_list_head();
                                neighbour_conductor != NULL;
                                    neighbour_conductor=neighbour_conductor->next)
                        {
                            // apply Newton iterators if we are in the same cell
                            if(current_cell == neighbour_cell)
                            {
                                if(current_newton_iterator < neighbour_newton_iterator)
                                {
                                    // get the neighbour particle
                                    Particle* neighbour_particle
                                            = neighbour_conductor->particle;

                                    compute_force(system_pt, molecule_pt,
                                                  current_particle,
                                                  neighbour_particle);
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
                                              current_particle,
                                              neighbour_particle);
                            }

                            // itreate the current conductor
                            // neighbour_conductor = neighbour_conductor->next;

                        } // end of loop over neighbour_conductor

                        // iterate the newton current newton iterator
                        current_newton_iterator++;

                        // itreate the current conductor
                        // current_conductor = current_conductor->next;

                    } // end of loop over current_conductor
                } // end of loop over n, number of neighbours
            } // end of loop over i, x-direction
        } // end of loop over j, y-direction
    } // end of loop over k, z-direction

    // set the reduction variables
    virial(0,0) = box_grad_00;
    virial(0,1) = box_grad_01;
    virial(1,1) = box_grad_11;

    if(number_of_cells_z != 0)
    {
        virial(0,2) = box_grad_02;
        virial(1,2) = box_grad_12;
        virial(2,2) = box_grad_22;
    }

    // printf("Force on particle 1: %f %f %f\n", molecule_pt->particle(0).f(0),molecule_pt->particle(0).f(1),molecule_pt->particle(0).f(2));
}
