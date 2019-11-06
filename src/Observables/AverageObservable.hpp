#ifndef AVERAGEOBSERVABLE_HPP
#define AVERAGEOBSERVABLE_HPP

using namespace::std;

/******************************************************************************
                              Average Observable
 *****************************************************************************/
class AverageObservable
{
public:
    AverageObservable();
    // destructor
    ~AverageObservable() {};
    // clear all the entries
    void clear();
    // make observation of the variable
    void observe(const double& value);
    // return the average of the observable
    double get_average();
    // get the instant last added value
    double get_instant();
    // return the variance of the observable
    double get_variance();
    // return the number of observations
    int get_nobservations();

private:
    double Average;
    double Average_sq;

    double n_observs;
    double instant_value;
};

#endif
