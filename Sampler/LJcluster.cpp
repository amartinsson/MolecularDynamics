#include <vector>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "BAOAB.hpp"
#include "Molecules.hpp"
#include "LennardJones.hpp"
#include "Array.hpp"

using namespace::std;
//
// double box_temp(std::vector<double>& Lp, double mass)
// {
//   unsigned dim = 3.0;
//
//   // mass and momentum for particle k
//   double momentum = 0.0;
//
//   // Temperature
//   double temperature = 0.0;
//
//   // loop over all the particles
//   for(unsigned j=0; j<dim; j++)
//   {
//       // dereference the momentum and mass
// 	  momentum = Lp[j];
// 	  temperature += momentum * momentum / mass;
//   }
//
//   // return the temperature
//   return 1.0 / dim * temperature;
// };

void print_positions(Molecule* molecule_pt, Matrix& S, const unsigned& time_stamp)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/Frames/frame_%i.csv", time_stamp);
  FILE* configuration_file = fopen(filename, "w");

  // dereference the dimension
  unsigned dim = molecule_pt->dim();
  unsigned number_of_particles = molecule_pt->nparticle();

  // make a particle pointer
//  Particle* particle_k = NULL;

  // print the correct stuff
  //for(unsigned k=0; k<number_of_particles; k++)
  // for(auto& part : molecule_pt->Particles)
  for(int i=0; i<number_of_particles; i++)
  {
    // dereference the particle
    Particle* part = molecule_pt->Particles.at(i);
    //particle_k = molecule_pt->particle_pt(k);

        // loop over the the number of dimensions
        for(unsigned j=0; j<dim; j++)
        {
            // print to file
            if(j == dim-1)
                fprintf(configuration_file, "%.4f", part->q(j));
            else
                fprintf(configuration_file, "%.4f, ", part->q(j));
        }

  	// dependedn on rigid body or not print different things
  	if(part->rigid_body())
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f, %.10f, %.10f, %.10f\n",
              part->Q(0,0), part->Q(0,1),
  		       part->Q(1,0), part->Q(1,1));
  	}
    else
    {
      // print end of line to file
      fprintf(configuration_file, "\n");
    }
  }
  // close the file
  fclose(configuration_file);

  // open the file to write to
  sprintf(filename, "Observables/fsimbox/frame_%i.csv", time_stamp);
  configuration_file = fopen(filename, "w");

  for(int i=0;i<10;i++)
  {
      if(i==0 | i==4)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", 0.0, 0.0, 0.0);
      else if(i==1)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,0), 0.0, 0.0);
      else if(i==2)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,0)+S(0,1), S(1,1), 0.0);
      else if(i==3)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,1), S(1,1), 0.0);
      else if(i==5 || i==9)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,2), S(1,2), S(2,2));
      else if(i==6)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,0)+S(0,2), S(1,2), S(2,2));
      else if(i==7)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,0)+S(0,1)+S(0,2), S(1,1)+S(1,2),S(2,2));
      else if(i==8)
        fprintf(configuration_file, "%.3f, %.3f, %.3f\n", S(0,1)+S(0,2), S(1,1)+S(1,2), S(2,2));
  }

  fclose(configuration_file);

};

void print_final_positions(Molecule* molecule_pt)
{

    // dereference the dimension
    unsigned dim = molecule_pt->dim();
    unsigned number_of_particles = molecule_pt->nparticle();

    // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/initial_position_%i.csv",
          number_of_particles);
  FILE* configuration_file = fopen(filename, "w");

  // make a particle pointer
  //Particle* particle_k = NULL;

  // print the correct stuff
  //for(unsigned k=0; k<number_of_particles; k++)
  // for(auto& part : molecule_pt->Particles)
  // {
  for(int i=0; i<number_of_particles; i++)
  {
        // dereference the particle
        Particle* part = molecule_pt->Particles.at(i);
    // dereference the particle
    //particle_k = molecule_pt->particle_pt(k);

    // loop over the the number of dimensions
    for(unsigned j=0; j<dim; j++)
    {
      // print to file
      if(j == dim-1)
        fprintf(configuration_file, "%.10f", part->q(j));
      else
        fprintf(configuration_file, "%.10f, ", part->q(j));
    }

  	// dependedn on rigid body or not print different things
  	if(part->rigid_body())
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f, %.10f, %.10f, %.10f\n",
              part->Q(0,0), part->Q(0,1),
  		        part->Q(1,0), part->Q(1,1));
  	}
    else
    {
      // print end of line to file
      fprintf(configuration_file, "\n");
    }
  }
  // close the file
  fclose(configuration_file);
};

void print_final_momentum(Molecule* molecule_pt)
{

    // dereference the dimension
    unsigned dim = molecule_pt->dim();
    unsigned number_of_particles = molecule_pt->nparticle();

    // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/initial_momentum_%i.csv",
          number_of_particles);
  FILE* configuration_file = fopen(filename, "w");

  // make a particle pointer
  //Particle* particle_k = NULL;

  // print the correct stuff
  //for(unsigned k=0; k<number_of_particles; k++)
  // for(auto& part : molecule_pt->Particles)
  // {
  for(int i=0; i<number_of_particles; i++)
  {
    // dereference the particle
    Particle* part = molecule_pt->Particles.at(i);
    // dereference the particle
    //particle_k = molecule_pt->particle_pt(k);

    // loop over the the number of dimensions
    for(unsigned j=0; j<dim; j++)
    {
      // print to file
      if(j == dim-1)
        fprintf(configuration_file, "%.10f", part->p(j));
      else
        fprintf(configuration_file, "%.10f, ", part->p(j));
    }

  	// dependedn on rigid body or not print different things
  	if(part->rigid_body())
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f\n", part->pi(0,0));
  	}
    else
    {
      // print end of line to file
      fprintf(configuration_file, "\n");
    }
  }
  // close the file
  fclose(configuration_file);
};

void print_temperature(const double& instant_temperature, double& temperature,
                       const double& time_stamp, const char* file_name,
                       const double& bmass, const double& blan)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%1.2f_%1.2f.csv", (name).c_str(),
            bmass, blan);
    FILE* temperature_file = fopen(filename, "a");

    fprintf(temperature_file, "%1.7e, %1.7e, %1.7e, %1.7e, %1.7e\n",
            time_stamp, bmass, blan, instant_temperature, temperature);

    // close the file
    fclose(temperature_file);
};

void print_box_temperature(const double& temperature, const double& time_stamp)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/box_temperature.csv");
  FILE* temperature_file = fopen(filename, "a");

  fprintf(temperature_file, "%1.4e, %1.3e\n", time_stamp, temperature);

    // close the file
  fclose(temperature_file);
};

void print_pressure(const double& instant_pressure, const double& pressure,
                    const double& time_stamp, const char* file_name,
                    const double& bmass, const double& blan)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%1.2f_%1.2f.csv", (name).c_str(),
            bmass, blan);
    FILE* pressure_file = fopen(filename, "a");

    fprintf(pressure_file, "%1.7e, %1.7e, %1.7e, %1.7e, %1.7e\n",
            time_stamp, bmass, blan, instant_pressure, pressure);

    // close the file
    fclose(pressure_file);
};

void print_density(const double& instant_volume, const unsigned& N,
                   const double& time_stamp, const unsigned& control_number)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/density_%d.csv", control_number);
  FILE* volume_file = fopen(filename, "a");

  //fprintf(volume_file, "%d, %.3f\n", time_stamp, volume);
  fprintf(volume_file, "%1.4e, %1.4e\n",
          time_stamp, (double)N / instant_volume);

    // close the file
  fclose(volume_file);
};


void print_volume(const double& volume, const double& time_stamp,
                  const double& bmass, const double& blan)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/volume_%1.2f_%1.2f.csv",
            bmass, blan);
  FILE* volume_file = fopen(filename, "a");

  fprintf(volume_file, "%1.5e, %1.5e, %1.5e, %1.5e\n", time_stamp,
            bmass, blan, volume);
  //fprintf(volume_file, "%1.4e, %.3f, %.3f, %.3f, %.3f\n", time_stamp, volume, box[0], box[1], box[2]);

    // close the file
  fclose(volume_file);
};

void print_final_box_size(const std::vector<double> box, const std::vector<double> momentum,
                            const unsigned& nparticles)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/box_size_%i.csv", nparticles);
  FILE* volume_file = fopen(filename, "w");

  //fprintf(volume_file, "%d, %.3f\n", time_stamp, volume);
  fprintf(volume_file, "%.7f, %.7f, %.7f, %.7f, %.7f, %.7f\n",
                       box[0], box[1], box[2],
                       momentum[0], momentum[1], momentum[2]);

    // close the file
  fclose(volume_file);
};

int main(int argc, char* argv[])
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    printf("Hello from process %d\n", world_rank);

    // parameters
    int SEED = 1234;
    unsigned dimension = 3;
    unsigned number_of_particles = 25*25;
    double temp = 0.2;

    double sigma = 1.0;
    double epsilon = 1.0;

    double sx = 30.0;
    double sy = 30.0;
    double sz = -1.0;

    double time_step = 0.001;
    double gamma = 1.0;
    double gamma_rot = 1.0;

    double target_pressure = 0.1;

    unsigned write_frequency = 100;

    double box_mass = 100.0;
    double box_mass_min = 10.0;
    double box_mass_max = 200.0;

    double npt_langevin_friction = 1.42;
    double npt_langevin_friction_min = 0.42;
    double npt_langevin_friction_max = 2.42;

    unsigned TIME = 1e2;
    bool rebuild_bool = false;

    unsigned npt_scheme_nr = 1;
    unsigned control_number = 0;

    double burn_in_fraction = 0.1;

  // cehck inouts
  for(int i=1;i<argc;i=i+2)
    {
        // Convert input to string
        std::string arg = argv[i];
        std::string arg2 = argv[i+1];
        if(arg == "--h")
        {
            double arg_in = std::stod(arg2);
            time_step = arg_in;
            printf("timestep set to: %.4f\n", time_step);
        }
        else if(arg == "--pressure")
        {
            double arg_in = std::stod(arg2);
            target_pressure = arg_in;
            printf("pressure set to: %1.2f\n", target_pressure);
        }
        else if(arg == "--temperature")
        {
            double arg_in = std::stod(arg2);
            temp = arg_in;
            printf("temperature set to: %1.2f\n", temp);
        }
        else if(arg == "--time")
        {
            unsigned arg_in = std::stod(arg2);
            TIME  = arg_in;
            printf("Total Simulation Time set to: %.0d\n", TIME);
        }
        else if(arg == "--nparticles")
        {
            unsigned arg_in = std::stod(arg2);
            number_of_particles  = arg_in;
            printf("number of particles set to: %.0d\n", number_of_particles);
        }
        else if(arg == "--S_x")
        {
            double arg_in = std::stod(arg2);
            sx = arg_in;
            printf("S_x set to: %1.2f\n", sx);
        }
        else if(arg == "--S_y")
        {
            double arg_in = std::stod(arg2);
            sy = arg_in;
            printf("S_y set to: %1.2f\n", sy);
        }
        else if(arg == "--S_z")
        {
            double arg_in = std::stod(arg2);
            sz = arg_in;
            printf("S_z set to: %1.2f\n", sz);
        }
        else if(arg == "--rbuild")
        {
            bool arg_in = std::stod(arg2);
            rebuild_bool  = arg_in;
            printf("Rebuilding.. %.0d\n", rebuild_bool);
        }
        else if(arg == "--bmass")
        {
            double arg_in = std::stod(arg2);
            box_mass = arg_in;
            printf("Thermal mass of box set to: %1.2f\n", box_mass);
        }
        else if(arg == "--bmass_max")
        {
            double arg_in = std::stod(arg2);
            box_mass_max = arg_in;
            printf("max Thermal mass of box set to: %1.2f\n", box_mass_max);
        }
        else if(arg == "--bmass_min")
        {
            double arg_in = std::stod(arg2);
            box_mass_min = arg_in;
            printf("min Thermal mass of box set to: %1.2f\n", box_mass_min);
        }
        else if(arg == "--blan")
        {
            double arg_in = std::stod(arg2);
            npt_langevin_friction = arg_in;
            printf("NPT friction of box set to: %1.2f\n",
                   npt_langevin_friction);
        }
        else if(arg == "--blan_max")
        {
            double arg_in = std::stod(arg2);
            npt_langevin_friction_max = arg_in;
            printf("max NPT friction of box set to: %1.2f\n",
                   npt_langevin_friction_max);
        }
        else if(arg == "--blan_min")
        {
            double arg_in = std::stod(arg2);
            npt_langevin_friction_min = arg_in;
            printf("min NPT friction of box set to: %1.2f\n",
                   npt_langevin_friction_min);
        }
        else if(arg == "--plan")
        {
            double arg_in = std::stod(arg2);
            gamma = arg_in;
            printf("Particle friction of simulation set to: %1.2f\n",
                   gamma);
        }
        else if(arg == "--burn")
        {
            double arg_in = std::stod(arg2);
            burn_in_fraction = arg_in;
            printf("burn in fraction set to: %1.2f\n", burn_in_fraction);
        }
        else if(arg == "--wfreq")
        {
            unsigned arg_in = std::stod(arg2);
            write_frequency = arg_in;
            printf("write frequency set to: %.0d\n", write_frequency);
        }
        else if(arg == "--nptnr")
        {
            unsigned arg_in = std::stod(arg2);
            npt_scheme_nr = arg_in;
            printf("NPT integration scheme set to: %.0d\n", npt_scheme_nr);
        }
        else if(arg == "--cntrl")
        {
            unsigned arg_in = std::stod(arg2);
            control_number = arg_in;
            printf("Control number set to: %.0d\n", control_number);
        }
        else if(arg == "--seed")
        {
            int arg_in = std::stod(arg2);
            SEED = arg_in;
            printf("Seed set to: %.0d\n", SEED);
        }

    }

    // initialise box size
    Matrix BoxS(2,2);
    if(sz > 0.0)
    {
        BoxS.resize(3,3);
        BoxS(2,2) = sz;
    }

    BoxS(0,0) = sx;
    BoxS(1,1) = sy;

    // pick the correct value for the box mass and langevin friction
    double dmass = (box_mass_max - box_mass_min) / sqrt((double)world_size);
    double dfric = (npt_langevin_friction_max - npt_langevin_friction_min)
                    / sqrt((double)world_size);

    double box_values[world_size];
    double lan_values[world_size];

    int sqrt_size = (int)sqrt((double)world_size);

    for(int i=0; i<sqrt_size; i++)
        for(int j=0; j<sqrt_size; j++)
        {
            box_values[j + i * sqrt_size] = box_mass_min + dmass * j;
            lan_values[j + i * sqrt_size] = npt_langevin_friction_min + dfric * i;
        }

    // add contorl number to seed differently
    SEED += control_number + (control_number + 1) * world_rank;

    box_mass = box_values[world_rank];
    npt_langevin_friction = lan_values[world_rank];

    printf("process %d has box_mass %1.3f\n",world_rank, box_mass);
    printf("process %d has npt_fric %1.3f\n",world_rank, npt_langevin_friction);

    // calculate the cut off
    double cut_off = 4.0 * pow(2.0, 1.0/6.0) * sigma;

    // calculate the number of steps
    unsigned number_of_steps = TIME / time_step;

    // calculate the number of burn in steps
    unsigned burn_in_steps = floor(burn_in_fraction * 1e3);

    // make particle standard values
    Vector q_0(dimension);
    Vector p_0(dimension);
    Matrix m_0(dimension, dimension);
    m_0.diag(1.0); // set mass to be one

    // make crystal molecule with number_of_particles
    Molecule* cluster = new Collection(q_0, p_0, m_0, temp,
                                       number_of_particles, dimension);

    // Make a Lennard Jones Force solver
    System* lennard_jones = new LennardJones(epsilon, sigma);

    // make integrator
    BAOAB* integrator = new BAOAB(1.0/temp, gamma, 0.0, time_step,
                                     lennard_jones, SEED);


    // set to evaluate with grid
    //integrator->integrate_with_grid(BoxS, cut_off, cluster);
    integrator->integrate_with_npt_grid(BoxS, cut_off, cluster,
                                        box_mass, target_pressure,
                                        npt_langevin_friction);

    for(unsigned i=0; i<number_of_steps; i++)
    {

        // integrate forward one step
        // double before = 0.0;
        // before = omp_get_wtime();

        integrator->integrate(cluster);
        if(integrator->cell_blow_up())
        {
            printf("ERROR: Breaking Process %d with mu = %1.3f, nu = %1.3f\n",
                    world_rank, box_mass, npt_langevin_friction);
            break;
        }

        // double after = omp_get_wtime();
        // printf("\tstep %4d took %1.4f\n", i, after-before);

        if(i % write_frequency == 0 && i != 0 && i > burn_in_steps)
        {
            double time_stamp = TIME * double(i) / double(number_of_steps);
            // print positions
            // Matrix simbox = integrator->npt_get_box();
            // print_positions(cluster, simbox, i / write_frequency);
            // update the pressure and temperature
            integrator->npt_update_pressure_temperature();
            // integrator->update_temperature();

            // print the pressure
            double instant_pressure = integrator->npt_get_instant_pressure();
            double pressure = integrator->npt_get_pressure();
            print_pressure(instant_pressure, pressure, time_stamp, "pressure",
                           box_mass, npt_langevin_friction);

            // printf("Pressure %.0d\t %1.5f\n", i, pressure);

            //print temperature
            double instant_temperature = integrator->npt_get_instant_temperature();
            double temperature = integrator->npt_get_temperature();
            print_temperature(instant_temperature, temperature, time_stamp,
                              "temperature",
                              box_mass, npt_langevin_friction);
            // printf("Temperature \t %1.5f\n", temperature);

            //print volume
            double volume = integrator->npt_get_volume();
            print_volume(volume, time_stamp,
                        box_mass,
                        npt_langevin_friction);
        }
    }

    // have mpi process wait if it blew up
    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize the MPI environment.
    MPI_Finalize();

} // end main
