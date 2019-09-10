#ifndef REPLICAEXCHANGEOBSERVABLE_HPP
#define REPLICAEXCHANGEOBSERVABLE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "ReplicaExchange.hpp"
#include "SystemTrajectory.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
template<class T>
class ReplicaExchangeObservable : public SystemObservable
{
public:
    ReplicaExchangeObservable(ReplicaExchange* method, T& ObsType,
        const int& recf, const int& rect);
    // empty destructor
    ~ReplicaExchangeObservable(){
        // delete obs;
    };
    // virtual function for updating
    void update();
    // get the average and instant vaalues
    double get_average() {/* empty */};
    // get the instant value
    double get_instant() {/* empty */};
    // print function at time_index
    void print(const char* file_name);

protected:

    ReplicaExchange* method;
    T* obs;
    unsigned N;
};

// constructor
template<class T>
ReplicaExchangeObservable<T>::ReplicaExchangeObservable(
    ReplicaExchange* method, T& ObsType,
        const int& recf, const int& rect)
            : N(method->replica.size()), SystemObservable(recf, rect),
                method(method)
{
    obs = (T*) malloc(N * sizeof(T));
};

template<class T>
void ReplicaExchangeObservable<T>::update()
{
    for(unsigned i=0; i<N; i++) {
        obs[i].update();
    }
}

// function for printing all the weights
template<class T>
void ReplicaExchangeObservable<T>::print(const char* file_name)
{
    for(unsigned i=0; i<N; i++) {
        obs[i].print(file_name, i);
    }
}

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
class RETrajectoryObservable : public SystemObservable
{
public:
    RETrajectoryObservable(ReplicaExchange* method, const int& recf,
        const int& rect);

    ~RETrajectoryObservable()
    {
        // for(unsigned i=0; i<N; i++)
        //     delete [] trajectory[i];
    };

    void update();
    // get the average and instant vaalues
    double get_average() {/* empty */};
    // get the instant value
    double get_instant() {/* empty */};
    // print the position
    void append_positions(const char* file_name, const double& time);
private:
    vector<SystemTrajectory*> trajectory;
    const unsigned N;
};

#endif
