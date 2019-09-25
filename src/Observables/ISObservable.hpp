#ifndef ISSOBSERVABLE_HPP
#define ISSOBSERVABLE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "InfiniteSwitch.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
class IsObservable : public SystemObservable
{
public:
    IsObservable(InfiniteSwitch* method, const int& recf,
                                 const int& rect);
    // empty destructor
    ~IsObservable(){};
    // virtual function for updating
    virtual void update() {/* empty */};
    // get the average and instant vaalues
    virtual double get_average() {/* empty */};
    // get the instant value
    virtual double get_instant() {/* empty */};

protected:
    InfiniteSwitch* method;
};

/******************************************************************************
                Simulated Tempering Temperature index Observable
 *****************************************************************************/
template<class T>
class IsWeightObservable : public IsObservable
{
public:
    // IsWeightObservable(InfiniteSwitchSimulatedTempering* method,
    //     SystemObservable* ObsType, const int& recf, const int& rect);
    IsWeightObservable(InfiniteSwitch* method,
        T& ObsType, const int& recf, const int& rect);
    // empty destructor
    ~IsWeightObservable() {};
    // virtual function for updating
    void update();
    // get the average and instant vaalues
    double get_average() {};
    // get the instant value
    double get_instant() {};
    // print function at time_index
    void print(const char* file_name);
    // holder for all the histograms
    // T* obs;
    vector<T> obs;

private:

    unsigned N;
};

// constructor
template<class T>
IsWeightObservable<T>::IsWeightObservable(
        InfiniteSwitch* method, T& ObsType,
            const int& recf, const int& rect)
                : IsObservable(method, recf, rect),
                    N(method->get_interpolation_points())
{
    // make break if the allocated memory isn't large enough,
    // dont know why this is.
    // cout << "size of ObsType: " << sizeof(ObsType) << endl;
    // obs = (T*) malloc(N * sizeof(ObsType));
    obs.reserve(N);
    // cout << "size of obs: " << sizeof(obs) << endl;
};

// function for updating
template<class T>
void IsWeightObservable<T>::update()
{
    // decode the observable weights
    vector<double> weight = method->get_observable_weights();

    for(unsigned i=0; i<N; i++) {
        obs[i].update(weight[i]);
    }
}

// function for printing all the weights
template<class T>
void IsWeightObservable<T>::print(const char* file_name)
{
    for(unsigned i=0; i<N; i++) {
        obs[i].print(file_name, i);
    }
}

#endif
