#include "SystemConfigurationalTemperature.hpp"

// constructor
SystemConfigurationalTemperature::SystemConfigurationalTemperature(
    Molecule* molecule_pt, const int& recf, const int& rect)
        : SystemTemperature(molecule_pt, recf, rect) {
            // keeper of average nabla square
            nablaSquare = new AverageObservable();

            // keeper of average laplace
            laplace = new AverageObservable();
        }

// destructor
SystemConfigurationalTemperature::~SystemConfigurationalTemperature() {
    // delete holders
    delete nablaSquare;
    delete laplace;
}

// get instant function
double SystemConfigurationalTemperature::get_instant()
{
    return nablaSquare->get_instant() / laplace->get_instant();
}

// get average function
double SystemConfigurationalTemperature::get_average()
{
    return nablaSquare->get_average() / laplace->get_average();
}

// update function
void SystemConfigurationalTemperature::update()
{
    if(recStep()) {
        // update the temperature
        update_temperature();
    }
}

// calculate the configurational temperature
void SystemConfigurationalTemperature::update_temperature() {
    // make temporary temperature
    double k = 0.0;
    
    for(const auto& particle : system->Particles) {
        k += particle.second->f.dot(particle.second->f);
    }

    // add to nabla square
    nablaSquare->observe(k);

    // add to the
    laplace->observe(SystemTemperature::system->laplace());
}
