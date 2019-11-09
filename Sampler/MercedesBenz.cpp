//
// Parameters that are stable:
// temperature = 0.05
// pressure = 0.1
// bmass =15
// S_x = 23
// nparticles = 256
// blan = 10.0




#include <vector>
#include <omp.h>
#include <stdio.h>
#include <math.h>

#include "BAOAB.hpp"
#include "Molecules.hpp"
#include "MercedesBenz.hpp"
#include "Array.hpp"
#include "SystemTrajectory.hpp"

#include "InfiniteSwitchSimulatedTempering.hpp"
#include "PressureInfiniteSwitchSimulation.hpp"

#include "SystemTrajectory.hpp"
#include "SystemTemperature.hpp"
#include "SystemEnergy.hpp"
#include "SystemPotentialEnergy.hpp"
#include "SystemConfigurationalTemperature.hpp"
#include "ISObservable.hpp"
#include "DISObservable.hpp"

using namespace::std;

int main(int argc, char* argv[])
{
    // Get the number of processes
    int world_size = 0;

    // Get the rank of the process
    int world_rank = 0;

    // parameters
    int SEED = 1234;
    unsigned dimension = 2;
    unsigned number_of_particles = 20*20;
    double temp = 0.1;

    double momIntertia = 0.01126;
    // double momIntertia = 0.001126;

    double sx = 30.0;
    double sy = 30.0;
    double sz = -1.0;

    double time_step = 0.001;
    double gamma = 1.0;
    double gamma_rot = 1.0;

    double target_pressure = 0.1;
    double pressure_max = 0.15;

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
        else if(arg == "--S_x")
        {
            double arg_in = std::stod(arg2);
            sx = arg_in;
            sy = sqrt(3.0) / 2.0 * sx;
            printf("S_x set to: %1.2f\n", sx);
            printf("S_y set to: %1.2f\n", sy);
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
        else if(arg == "--rlan")
        {
            double arg_in = std::stod(arg2);
            gamma_rot = arg_in;
            printf("Particle rotation friction of simulation set to: %1.2f\n",
                   gamma_rot);
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
        else if(arg == "--seed")
        {
            int arg_in = std::stod(arg2);
            SEED = arg_in;
            printf("Seed set to: %.0d\n", SEED);
        }
        else if(arg == "--rbuild")
        {
            bool arg_in = std::stod(arg2);
            rebuild_bool  = true;
            printf("Rebuilding.. %.0d\n", rebuild_bool);
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

    // add contorl number to seed differently
    SEED += control_number + (control_number + 1) * world_rank;

    // calculate the number of steps
    unsigned number_of_steps = TIME / time_step;

    // calculate the number of burn in steps
    unsigned burn_in_steps = floor(burn_in_fraction * number_of_steps);

    // make particle standard values
    Vector q_0(dimension);
    Vector p_0(dimension);
    Matrix Q_0(dimension, dimension);
    Matrix pi_0(1, 1);
    Matrix m_0(dimension, dimension);
    Matrix I_0(dimension, dimension);

    // set mass and moment of inertia
    m_0.diag(1.0);
    I_0.diag(momIntertia);

    // make crystal molecule with number_of_particles
    unsigned narms = 3;
    Molecule* cluster = new AniCollection(q_0, p_0, Q_0, pi_0, m_0, I_0, narms,
                                          temp, number_of_particles, dimension);

    // Make a Mercedes Benz Force solver
    double sigma = 0.7;
    double epsilon_hb = -1.0;
    double sigma_hb = 0.085;
    double r_hb = 1.0;
    double epsilon = 0.1 * fabs(epsilon_hb);
    // double epsilon = -0.1;
    // double epsilon = 0.1; // used in the MB paper
    // double epsilon = 0.4; // used in the best paper

    System* mercedes_benz = new MercedesBenz(epsilon, sigma, epsilon_hb, sigma_hb, r_hb);

    // make integrator
    BAOAB* integrator = new BAOAB(1.0/temp, gamma, gamma_rot, time_step, mercedes_benz, SEED);


    // calculate the cut off
    double cut_off = 4.0 * r_hb;

    // set the boundaries for the rmin rmax calculations
    double rmin = 0.0;
    double rmax = cut_off;
    int nrdist = 200;

    // integrate with npt grid
    integrator->integrate_with_npt_grid(BoxS, cut_off, cluster, box_mass, target_pressure, npt_langevin_friction, write_frequency, burn_in_steps);

printf("about to read\n");
    if(rebuild_bool == true)
    {
        char position[50];
        char volume[50];
        char momentum[50];

        sprintf(position, "Observables/Initial/position_0.csv");
        sprintf(volume, "Observables/Initial/volume_0.csv");
        sprintf(momentum, "dummy");
        integrator->npt_set_initial(cluster, position, momentum, volume);

        printf("Read in:\n\t%2.3f,%2.3f\n\t%2.3f,%2.3f\n",
        integrator->npt_obj().S(0,0), integrator->npt_obj().S(0,1), integrator->npt_obj().S(1,0), integrator->npt_obj().S(1,1));
    }

    // calculate the order parameter
    integrator->npt_obj().set_to_calculate_sphere_order_param(cluster, 6, 1.225 * sigma);

    // calculate the radial distribution
    integrator->npt_obj().set_to_calculate_radial_dist(rmin, rmax, nrdist);

    // make system trajectory object
    SystemTrajectory traj = SystemTrajectory(cluster);

    traj.set_simbox(integrator->npt_obj().S);

    // energy recorder
    SystemEnergy energy = SystemEnergy(cluster, write_frequency, burn_in_steps);

    // set the particles on a hexagonal grid
    integrator->npt_obj().set_particles_hexagonal(cluster);

    traj.print_positions("frame", 0);
    traj.print_simbox("simbox", 0);
    // exit(-1);

    for(unsigned i=0; i<number_of_steps + burn_in_steps; i++)
    {
        // integrate forward
        integrator->integrate(cluster);


        // printf("V = %f\n", cluster->potential());
        // update pressure temperature// update pressure temperature
        integrator->npt_obj().update_pressure_temperature();

        if(integrator->cell_blow_up())
        {
            printf("step %d\n", i);
            printf("ERROR: Breaking Process %d with mu = %1.3f, nu = %1.3f\n",
                    world_rank, box_mass, npt_langevin_friction);
            break;
        }


        double tForceX = 0.0;
        double tForceY = 0.0;

        double tTau = 0.0;

        for(auto& particle : cluster->Particles) {
            tForceX += particle.second->f(0);
            tForceX += particle.second->f(1);

            tTau += particle.second->tau(0,0);
        }

        // printf("Fx = %1.3f, Fy = %1.3f, Tau=%1.3f\n", tForceX, tForceY, tTau);
        //
        // printf("Particle 0 = [%1.3e, %1.3e]\n", cluster->particle(0).f(0), cluster->particle(0).f(1));

        // update the energy
        energy.update_potential();
        if(i % 1000 * write_frequency == 0) {

            traj.print_positions("frame", i/write_frequency);
            traj.print_simbox("simbox", i/write_frequency);
        }

        // exit(-1);

        // if(i % write_frequency == 0 && i != 0 && i > burn_in_steps)
        if(i % write_frequency == 0 && i > burn_in_steps)
        {
            double time_stamp = TIME * double(i - burn_in_steps) / double(number_of_steps);


            // print temperature
            integrator->npt_obj().Temperature_pt->print("temperature",
                                                time_stamp,
                                                control_number + world_rank);
            // print pressure
            integrator->npt_obj().Pressure_pt->print("pressure",
                                                time_stamp,
                                                control_number + world_rank);
            // print volume
            integrator->npt_obj().Volume_pt->print("volume", time_stamp,
                                                control_number + world_rank);
            // print volume
            integrator->npt_obj().SphereOrder_pt->print("order", time_stamp,
                                                control_number + world_rank);
            // print potential energy
            energy.print("potential", time_stamp, control_number + world_rank);
        }
    }

    // print the distribution function
    integrator->npt_obj().Radial_pt->print("rdist", 0.0,
                                                control_number + world_rank);

} // end main
