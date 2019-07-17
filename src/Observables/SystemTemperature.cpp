#include "SystemTemperature.hpp"

// constructor
SystemTemperature::SystemTemperature(Molecule* molecule_pt, const int& recf,
                                     const int& rect)
    : SystemObservable(recf, rect), system(molecule_pt)
{
    // average holder for the temperature
    Temp = new AverageObservable();
}

// get instant function
double SystemTemperature::get_instant()
{
    return Temp->get_instant();
}

// get instant function
double SystemTemperature::get_average()
{
    return Temp->get_average();
}

// update function
void SystemTemperature::update()
{
    if(recStep()) {
        // update the temperature
        update_temperature();
    }
}

// calculate momentum temperature
void SystemTemperature::update_temperature()
{
    // make an average observable
    AverageObservable temp = AverageObservable();

    // loop overa all the particles
    for(const auto& particle : system->Particles)
    {
        temp.observe(particle.second->p.dot(
                     particle.second->m.inv() * particle.second->p));

        // if(particle.second->rigid_body())
        //     temp.observe(particle.second->pi(0, 0) *
        //                  particle.second->I(0, 0) * particle.second->pi(0,0));

    }

    // add the observed temperature to the mTemp observable
    double dim = double(system->dim());

    // if(system->particle(0).rigid_body())
    //     dim++;

    double N = double(system->nparticle());

    if(N > 1.0)
        Temp->observe(N * temp.get_average() / (dim * (N - 1.0)));
    else {
        Temp->observe(temp.get_average() / dim);
    }
}
