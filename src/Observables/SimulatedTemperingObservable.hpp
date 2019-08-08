#ifndef SIMULATEDTEMPERINGOBSERVABLE_HPP
#define SIMULATEDTEMPERINGOBSERVABLE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "SimulatedTempering.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
class SimulatedTemperingObservable : public SystemObservable
{
public:
    SimulatedTemperingObservable(SimulatedTempering* method, const int& recf,
                                 const int& rect);
    // empty destructor
    ~SimulatedTemperingObservable(){};
    // virtual function for updating
    virtual void update() {/* empty */};
    // get the average and instant vaalues
    virtual double get_average() {/* empty */};
    // get the instant value
    virtual double get_instant() {/* empty */};

protected:
    SimulatedTempering* method;
};

/******************************************************************************
                Simulated Tempering Temperature index Observable
 *****************************************************************************/
class SimTempTemperatureObservable : public SimulatedTemperingObservable
{
public:
    SimTempTemperatureObservable(SimulatedTempering* method);
    // empty destructor
    ~SimTempTemperatureObservable(){};
    // virtual function for updating
    void update();
    // get the average and instant vaalues
    double get_average();
    // get the instant value
    double get_instant();
    // print function at time_index
    void print(const char* file_name, const double& time, const unsigned& index);

private:
    AverageObservable temp_target;
    AverageObservable temp_index;
};

#endif
