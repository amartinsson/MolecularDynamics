#ifndef DISSOBSERVABLE_HPP
#define DISSOBSERVABLE_HPP

#include "AverageObservable.hpp"
#include "SystemObservable.hpp"
#include "DoubleInfiniteSwitch.hpp"
#include "Array.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
class DisObservable : public SystemObservable
{
public:
    DisObservable(DoubleInfiniteSwitch* method, const int& recf,
                                 const int& rect);
    // empty destructor
    ~DisObservable(){};
    // virtual function for updating
    virtual void update() {/* empty */};
    // get the average and instant vaalues
    virtual double get_average() {/* empty */};
    // get the instant value
    virtual double get_instant() {/* empty */};

protected:
    DoubleInfiniteSwitch* method;
};

/******************************************************************************
                Simulated Tempering Temperature index Observable
 *****************************************************************************/
template<class T>
class DisWeightObservable : public DisObservable
{
public:
    // DisWeightObservable(InfiniteSwitchSimulatedTempering* method,
    //     SystemObservable* ObsType, const int& recf, const int& rect);
    DisWeightObservable(DoubleInfiniteSwitch* method,
        T& ObsType, const int& recf, const int& rect);
    // empty destructor
    ~DisWeightObservable() {};
    // virtual function for updating
    void update();
    // get the average and instant vaalues
    double get_average() {};
    // get the instant value
    double get_instant() {};
    // print function at time_index
    void print(const char* file_name);
    // holder for all the histograms
    vector<T*> obs;

private:

    unsigned None;
    unsigned Ntwo;
};

// constructor
template<class T>
DisWeightObservable<T>::DisWeightObservable(
        DoubleInfiniteSwitch* method, T& ObsType,
            const int& recf, const int& rect)
                : DisObservable(method, recf, rect),
                    None(method->get_interpolation_points()[0]),
                        Ntwo(method->get_interpolation_points()[1])
{
    for(unsigned i=0; i<None; i++) {
        obs.push_back((T*) malloc(Ntwo * sizeof(T)));
    }
};

// function for updating
template<class T>
void DisWeightObservable<T>::update()
{
    // decode the observable weights
    Matrix weight = method->get_observable_weights();

#pragma omp for collapse(2)
    for(unsigned i=0; i<None; i++) {
        for(unsigned j=0; j<Ntwo; j++) {
            obs[i][j].update(weight(i,j));
        }
    }
}

// function for printing all the weights
template<class T>
void DisWeightObservable<T>::print(const char* file_name)
{

#pragma omp for collapse(2)
    for(unsigned i=0; i<None; i++) {
        for(unsigned j=0; j<Ntwo; j++) {
            char filename[50];

            string name(file_name);

            sprintf(filename, "%s_%i", (name).c_str(), i);
            
            obs[i][j].print(filename, j);
        }
    }
}

#endif
