#include <vector>
#include <omp.h>
#include <stdio.h>
#include <math.h>

#include "BAOAB.hpp"
#include "Molecules.hpp"
#include "CurieWeiss.hpp"
#include "InfiniteSwitchSimulatedTempering.hpp"
#include "SimulatedTempering.hpp"

/* Observables to include */
#include "SystemTrajectory.hpp"
#include "SystemTemperature.hpp"
#include "SystemEnergy.hpp"
#include "SystemPotentialEnergy.hpp"
#include "SystemConfigurationalTemperature.hpp"
#include "MagnetisationObservable.hpp"

int main(int argc, char* argv[])
{
    // parameters
    int SEED = 1234;
    unsigned dimension = 1;
    double temp = 0.25;

    double b = 1.0;
    unsigned number_of_particles = 10;

    double time_step = 0.1;
    double gamma = 1.0;
    unsigned write_frequency = 1;

    unsigned TIME = 1e2;
    double burn_in_fraction = 0.0;

    double px = 0.0;
    double qx = 0.0;

    double temp_max = 2.0;

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
        else if(arg == "--temperature")
        {
            double arg_in = std::stod(arg2);
            temp = arg_in;
            printf("temperature set to: %1.2f\n", temp);
        }
        else if(arg == "--temp_max")
        {
            double arg_in = std::stod(arg2);
            temp_max = arg_in;
            printf("max temperature set to: %1.2f\n", temp_max);
        }
        else if(arg == "--time")
        {
            unsigned arg_in = std::stod(arg2);
            TIME  = arg_in;
            printf("Total Simulation Time set to: %.0d\n", TIME);
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
        else if(arg == "--seed")
        {
            int arg_in = std::stod(arg2);
            SEED = arg_in;
            printf("Seed set to: %.0d\n", SEED);
        }
        else if(arg == "--npart")
        {
            int arg_in = std::stod(arg2);
            number_of_particles = arg_in;
            printf("Number of particles set to: %.0d\n", number_of_particles);
        }
        else if(arg == "--px")
        {
            double arg_in = std::stod(arg2);
            px = arg_in;
        }
        else if(arg == "--qx")
        {
            double arg_in = std::stod(arg2);
            qx = arg_in;
            printf("Starting at: %1.2f\n", qx);
        }
        else if(arg == "--b")
        {
            double arg_in = std::stod(arg2);
            b = arg_in;
            printf("external field set to: %1.4f\n", b);
        }
    }

    // calculate the number of steps
    unsigned number_of_steps = TIME / time_step;

    // calculate the number of burn in steps
    unsigned burn_in_steps = floor(burn_in_fraction * number_of_steps);

    // make particle standard values
    Vector q_0(dimension);
    Vector p_0(dimension);
    Matrix m_0(dimension, dimension);

    // set mass and position
    // m_0.diag(1.0);
    m_0(0,0) = 1.0;
    // q_0(0) = -1.0;
    // q_0(1) = -1.0;
    q_0(0) = qx;

    p_0(0) = px;

    // make crystal molecule with number_of_particles
    Molecule* cluster = new Collection(q_0, p_0, m_0, temp, number_of_particles, dimension);

    // Make a Lennard Jones Force solver
    System* curie_weiss = new CurieWeiss(b, number_of_particles);
    // solve for the force
    curie_weiss->compute_force(cluster);

    // make integrator
    BAOAB* integrator = new BAOAB(1.0/temp, gamma, 0.0, time_step,
                                  curie_weiss, SEED);

    // make system trajectory object
    SystemTrajectory traj =
        SystemTrajectory(cluster);

    SystemConfigurationalTemperature stemp =
        SystemConfigurationalTemperature(cluster, write_frequency,
                                         burn_in_steps);
    SystemEnergy energy =
        SystemEnergy(cluster, write_frequency, burn_in_steps);

    MagnetisationObservable magnet =
        MagnetisationObservable(cluster, -1.0, 1.0, 100,
                                    write_frequency, burn_in_steps);

    // set to integrate with isst
    unsigned nint = 10;
    double isstTau = 1.0;
    integrator->integrate_with_isst(cluster, temp, temp_max, nint, time_step, isstTau);

    for(unsigned i=0; i<number_of_steps; i++) {
        // integrate forward
        integrator->integrate(cluster);

        // append to trajectory
        double t = (double)i / (double)number_of_steps * TIME;
        traj.append_positions("trajectory", t);

        // update the temperature
        stemp.update();
        stemp.print("temperature", t);

        // update the energy
        energy.update();
        energy.print("energy", t);

        // update the magnetisation
        magnet.update();
    }

    // print the magnetisation
    magnet.print("magnetisation", 0);
} // end main
