#include "SimulatedTemperingObservable.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
// constructor for abstract base class
SimulatedTemperingObservable::SimulatedTemperingObservable(
    SimulatedTempering* method, const int& recf, const int& rect)
        : SystemObservable(recf, rect), method(method) {};


/******************************************************************************
                Simulated Tempering Temperature index Observable
 *****************************************************************************/
SimTempTemperatureObservable::SimTempTemperatureObservable(
    SimulatedTempering* method)
        : SimulatedTemperingObservable(method, 1, 0), temp_target(),
            temp_index() {}

// virtual function for updating
void SimTempTemperatureObservable::update() {
    // update on every step
    double index = (double)method->get_temperature_index();
    double temp = method->get_temperature();

    // update the averages'
    temp_index.observe(index);
    temp_target.observe(temp);
}

// get the average and instant vaalues
double SimTempTemperatureObservable::get_average() {
    return temp_target.get_average();
}

// get the instant value
double SimTempTemperatureObservable::get_instant() {
    return temp_target.get_instant();
}

// print function at time_index
void SimTempTemperatureObservable::print(const char* file_name,
                                         const double& time,
                                         const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    // print the current state and long term average
    fprintf(file, "%1.7e, %.0d, %1.3e, %1.3e\n", time,
            (int)temp_index.get_instant(), temp_target.get_instant(),
            temp_target.get_average());

    // close the file
    fclose(file);
}
