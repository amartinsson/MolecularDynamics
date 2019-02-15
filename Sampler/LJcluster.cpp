#include <vector>
#include <omp.h>
#include <stdio.h>
#include <math.h>
// #include <mpi.h>

#include "BAOAB.hpp"
#include "Molecules.hpp"
#include "LennardJones.hpp"

using namespace::std;

double box_temp(std::vector<double>& Lp, double mass)
{
  unsigned dim = 3.0;

  // mass and momentum for particle k
  double momentum = 0.0;

  // Temperature
  double temperature = 0.0;

  // loop over all the particles
  for(unsigned j=0; j<dim; j++)
  {
      // dereference the momentum and mass
	  momentum = Lp[j];
	  temperature += momentum * momentum / mass;
  }

  // return the temperature
  return 1.0 / dim * temperature;
};

void print_positions(Molecule* molecule_pt, const unsigned& time_stamp)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/Frames/frame_%i.csv", time_stamp);
  FILE* configuration_file = fopen(filename, "w");

  // dereference the dimension
  unsigned dim = molecule_pt->dim();
  unsigned number_of_particles = molecule_pt->nparticle();

  // make a particle pointer
  Particle* particle_k = NULL;

  // print the correct stuff
  for(unsigned k=0; k<number_of_particles; k++)
  {
    // dereference the particle
    particle_k = molecule_pt->particle_pt(k);

    // loop over the the number of dimensions
    for(unsigned j=0; j<dim; j++)
    {
      // print to file
      if(j == dim-1)
        fprintf(configuration_file, "%.4f", *particle_k->q_pt(j));
      else
        fprintf(configuration_file, "%.4f, ", *particle_k->q_pt(j));
    }

  	// dependedn on rigid body or not print different things
  	if(particle_k->rigid_body() != false)
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f, %.10f, %.10f, %.10f\n",
              *particle_k->Q_pt(0,0), *particle_k->Q_pt(0,1),
  		        *particle_k->Q_pt(1,0), *particle_k->Q_pt(1,1));
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
  Particle* particle_k = NULL;

  // print the correct stuff
  for(unsigned k=0; k<number_of_particles; k++)
  {
    // dereference the particle
    particle_k = molecule_pt->particle_pt(k);

    // loop over the the number of dimensions
    for(unsigned j=0; j<dim; j++)
    {
      // print to file
      if(j == dim-1)
        fprintf(configuration_file, "%.10f", *particle_k->q_pt(j));
      else
        fprintf(configuration_file, "%.10f, ", *particle_k->q_pt(j));
    }

  	// dependedn on rigid body or not print different things
  	if(particle_k->rigid_body() != false)
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f, %.10f, %.10f, %.10f\n",
              *particle_k->Q_pt(0,0), *particle_k->Q_pt(0,1),
  		        *particle_k->Q_pt(1,0), *particle_k->Q_pt(1,1));
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
  Particle* particle_k = NULL;

  // print the correct stuff
  for(unsigned k=0; k<number_of_particles; k++)
  {
    // dereference the particle
    particle_k = molecule_pt->particle_pt(k);

    // loop over the the number of dimensions
    for(unsigned j=0; j<dim; j++)
    {
      // print to file
      if(j == dim-1)
        fprintf(configuration_file, "%.10f", *particle_k->p_pt(j));
      else
        fprintf(configuration_file, "%.10f, ", *particle_k->p_pt(j));
    }

  	// dependedn on rigid body or not print different things
  	if(particle_k->rigid_body() != false)
  	{
      // print the rotation matrix to the file
      fprintf(configuration_file, ", %.10f\n", *particle_k->pi_pt(0));
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
                       const unsigned& control_number)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), control_number);
    FILE* temperature_file = fopen(filename, "a");

    fprintf(temperature_file, "%1.7e, %1.7e, %1.7e\n",
            time_stamp, instant_temperature, temperature);

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
                    const unsigned& control_number)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), control_number);
    FILE* pressure_file = fopen(filename, "a");

    fprintf(pressure_file, "%1.7e, %1.7e, %1.7e\n",
            time_stamp, instant_pressure, pressure);

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


void print_volume(const double& volume, const std::vector<double> box, const double& time_stamp)
{
  // open the file to write to
  char filename[50];
  sprintf(filename, "Observables/volume.csv");
  FILE* volume_file = fopen(filename, "a");

  //fprintf(volume_file, "%d, %.3f\n", time_stamp, volume);
  fprintf(volume_file, "%1.4e, %.3f, %.3f, %.3f, %.3f\n", time_stamp, volume, box[0], box[1], box[2]);

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
    //MPI_Init(NULL, NULL);

    // Get the number of processes
    //int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank = 0;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    // char processor_name[MPI_MAX_PROCESSOR_NAME];
    // int name_len;
    // MPI_Get_processor_name(processor_name, &name_len);
    //
    // // Print off a hello world message
    // printf("Hello world from processor %s, rank %d out of %d processors\n",
    //        processor_name, world_rank, world_size);

    // print GPU information
    // printf("Number of devices detected: %.d\n",  omp_get_num_devices());
    // omp_set_default_device(1);
    // printf("The Default device is: %.d\n", omp_get_default_device());

//     std::cout << omp_get_num_devices() << " number of devices\n";
  // parameters
  int SEED = 1234;
  unsigned dimension = 2;
  unsigned number_of_particles = 25*25;
  double temp = 0.2;

  double sigma = 0.5;
  double epsilon = 0.01;

  double a_x = 30.0;
  double b_x = 0.0;
  double b_y = 30.0;

  double time_step = 0.001;
  double gamma = 1.0;
  double gamma_rot = 1.0;

  double target_pressure = 0.1;

  unsigned write_frequency = 100;

  double box_mass = 100.0;
  double npt_langevin_friction = 1.42;

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
        else if(arg == "--a_x")
        {
            double arg_in = std::stod(arg2);
            a_x = arg_in;
            printf("a_x set to: %1.2f\n", a_x);
        }
        else if(arg == "--b_x")
        {
            double arg_in = std::stod(arg2);
            b_x = arg_in;
            printf("b_x set to: %1.2f\n", b_x);
        }
        else if(arg == "--b_y")
        {
            double arg_in = std::stod(arg2);
            b_y = arg_in;
            printf("b_y set to: %1.2f\n", b_y);
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
        else if(arg == "--blan")
        {
            double arg_in = std::stod(arg2);
            npt_langevin_friction = arg_in;
            printf("NPT friction of box set to: %1.2f\n",
                   npt_langevin_friction);
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

    // add contorl number to seed differently
    SEED += control_number + (control_number + 1) * world_rank;

    // calculate the cut off
    double cut_off = 4.0 * pow(2.0, 1.0/6.0) * sigma;

    // calculate the number of steps
    unsigned number_of_steps = TIME / time_step;

    // calculate the number of burn in steps
    unsigned burn_in_steps = floor(burn_in_fraction * number_of_steps);

    // make particle standard values
    vector<double> q_0(dimension, 0.0); // Initial position
    vector<double> p_0(dimension, 0.0); // Initial momentum
    vector<double> f_0(dimension, 0.0); // Initial Force
    vector<double> M(dimension, 1.0);   // Mass

    // make crystal molecule with number_of_particles
    Molecule* cluster = new Crystal(q_0, p_0, M, number_of_particles, temp);

    // Make a Lennard Jones Force solver
    System* lennard_jones = new LennardJones(epsilon, sigma);

    // make integrator
    BAOAB* integrator = new BAOAB(1.0/temp, gamma, 0.0, time_step,
                                     lennard_jones, SEED);

    // set to evaluate with grid
    integrator->integrate_with_npt_grid(a_x, b_x, b_y, cut_off, cluster,
                                        box_mass, target_pressure,
                                        npt_langevin_friction);

    // set the integrator scheme
    integrator->set_npt_integrator_version(npt_scheme_nr);

    if(rebuild_bool)
        integrator->npt_set_initial(cluster, "Init/initial_position.csv",
                                    "Init/initial_momentum.csv",
                                    "Init/initial_volume.csv");

    // char filename[50];
    // sprintf(filename, "timing_0.csv");
    // FILE* output = fopen(filename, "w");

    // Integration loop
// #pragma omp parallel shared(integrator, cluster) firstprivate(lennard_jones)
{
    for(unsigned i=0; i<number_of_steps; i++)
    {
        // integrate forward one step
        double before = 0.0;
        // #pragma omp master
        before = omp_get_wtime();

        integrator->integrate(cluster);

        // #pragma omp master
        {
            double after = omp_get_wtime();
            // fprintf(output, "%f\n", after-before);
            printf("%f\n", after-before);
        }



        // print the position
// #pragma omp single
        {
            if(i % write_frequency == 0 && i != 0 && i > burn_in_steps)
            {
                double time_stamp = TIME * double(i) / double(number_of_steps);

                // update the pressure and temperature
                integrator->npt_update_pressure_temperature();
                // integrator->update_temperature();

                // print the pressure
                double instant_pressure = integrator->npt_get_instant_pressure();
                double pressure = integrator->npt_get_pressure();
                print_pressure(instant_pressure, pressure, time_stamp, "pressure",
                               control_number + (control_number + 1) * world_rank);

                // printf("Pressure on step %.0d is %1.3f\n", i, pressure);

                // print temperature
                double instant_temperature = integrator->npt_get_instant_temperature();
                double temperature = integrator->npt_get_temperature();
                print_temperature(instant_temperature, temperature, time_stamp,
                                  "temperature", control_number + (control_number + 1) * world_rank);

                // //print the density
                // double volume = integrator->get_instant_npt_volume();
                // print_density(volume, number_of_particles, time_stamp, control_number);
                //
                // //double volume = integrator->get_npt_volume();
                // std::vector<double> box_coordinate = integrator->get_npt_box_coord();
                // print_volume(volume, box_coordinate, time_stamp);

                // print the positions
                // print_positions(cluster, i / write_frequency);
            }

            // only print on the last step
            // if(i == number_of_steps-1)
            // {
            //     // get box size and momentum
            //     std::vector<double> box_coordinate = integrator->get_npt_box_coord();
            //     std::vector<double> box_momentum = integrator->get_npt_box_momentum();
            //     print_final_box_size(box_coordinate, box_momentum, number_of_particles);
            //
            //     // print the position of all the particles
            //     print_final_positions(cluster);
            //
            //     // print the momentum of all the particles
            //     print_final_momentum(cluster);
            // }
        }
    }
}

    // Finalize the MPI environment.
    // MPI_Finalize();

} // end main
