#include "AverageObservable.hpp"

using namespace::std;

/******************************************************************************
                              Average Observable
 *****************************************************************************/
AverageObservable::AverageObservable()
{
    Average = 0.0;
    Average_sq = 0.0;
    n_observs = 0.0;
    instant_value = 0.0;
}

// destructor
AverageObservable::~AverageObservable()
{
    delete &Average;
    delete &Average_sq;
    delete &n_observs;
    delete &instant_value;
}

// make observation
void AverageObservable::observe(const double& value)
{
    #pragma omp critical
    {
        // bump up the observation counter
        n_observs++;
        double N = n_observs;

        // sum up the values
        Average = value / N + (N - 1.0) / N * Average;
        // sum up the variance
        Average_sq = value * value / N + (N-1.0) / N * Average_sq;
        // update the cached instant value
        instant_value = value;
    }
};

// get the average of the observable
double AverageObservable::get_average()
{
    // return the average
    return Average;
};

// get the instant last added value
double AverageObservable::get_instant()
{
    // return the instant value
    return instant_value;
}

// get the variance of the observable
double AverageObservable::get_variance()
{
    // return the variance
    return Average_sq - Average * Average;
};
