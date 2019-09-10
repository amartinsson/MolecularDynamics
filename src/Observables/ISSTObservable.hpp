#ifndef ISSTOBSERVABLE_HPP
#define ISSTOBSERVABLE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "InfiniteSwitchSimulatedTempering.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
class IsstObservable : public SystemObservable
{
public:
    IsstObservable(InfiniteSwitchSimulatedTempering* method, const int& recf,
                                 const int& rect);
    // empty destructor
    ~IsstObservable(){};
    // virtual function for updating
    virtual void update() {/* empty */};
    // get the average and instant vaalues
    virtual double get_average() {/* empty */};
    // get the instant value
    virtual double get_instant() {/* empty */};

protected:
    InfiniteSwitchSimulatedTempering* method;
};

/******************************************************************************
                Simulated Tempering Temperature index Observable
 *****************************************************************************/
template<class T>
class IsstWeightObservable : public IsstObservable
{
public:
    // IsstWeightObservable(InfiniteSwitchSimulatedTempering* method,
    //     SystemObservable* ObsType, const int& recf, const int& rect);
    IsstWeightObservable(InfiniteSwitchSimulatedTempering* method,
        T& ObsType, const int& recf, const int& rect);
    // empty destructor
    ~IsstWeightObservable() {};
    // virtual function for updating
    void update();
    // get the average and instant vaalues
    double get_average() {};
    // get the instant value
    double get_instant() {};
    // print function at time_index
    void print(const char* file_name);
    // holder for all the histograms
    T* obs;

private:

    unsigned N;
};

// constructor
template<class T>
IsstWeightObservable<T>::IsstWeightObservable(
        InfiniteSwitchSimulatedTempering* method, T& ObsType,
            const int& recf, const int& rect)
                : IsstObservable(method, recf, rect),
                    N(method->get_interpolation_points())
{
    obs = (T*) malloc(N * sizeof(T));
};

// function for updating
template<class T>
void IsstWeightObservable<T>::update()
{
    // decode the observable weights
    vector<double> weight = method->get_observable_weights();

    for(unsigned i=0; i<N; i++) {
        obs[i].update(weight[i]);
    }
}

// function for printing all the weights
template<class T>
void IsstWeightObservable<T>::print(const char* file_name)
{
    for(unsigned i=0; i<N; i++) {
        obs[i].print(file_name, i);
    }
}

#endif
