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


// calculate the configurational temperature
void SystemConfigurationalTemperature::update_temp() {

}
