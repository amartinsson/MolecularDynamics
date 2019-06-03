#ifndef HISTOBSERVABLE_HPP
#define HISTOBSERVABLE_HPP

#include <gsl/gsl_histogram.h>
#include <omp.h>

using namespace::std;

/******************************************************************************
                        Wrapper for the gsl histogram observable
 *****************************************************************************/
class HistObservable
{
public:
    // constructor
    HistObservable(const double& low, const double& up, const int& n);
    // empty destructor
    ~HistObservable()
    {
        // free the memory for the histogram
        gsl_histogram_free(histogram);
    };
    // virtual function for updating
    void observe(const double& value);
    void observe(const double& value, const double& weight);
    // get the vaolue of bin i
    double get_value(const int& i);
    // get the number of bins
    int get_nbins();
    // get the center of the ith bin
    double get_bin_center(const int& i);
private:
    // gsl histogram structure
    gsl_histogram* histogram;
    // number of bins
    int N;
};
#endif
