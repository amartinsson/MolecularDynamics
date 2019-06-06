#include "HistObservable.hpp"

using namespace::std;

// constructor
HistObservable::HistObservable(const double& low, const double& up,
                               const int& n) : N(n)
{
    // allocate memory for histogram
    histogram = gsl_histogram_alloc(n);
    // make the bins uniform between low and up
    gsl_histogram_set_ranges_uniform(histogram, low, up);
}

// add an observation to the histogram
void HistObservable::observe(const double& value)
{
    // increment histogram with one
    #pragma omp critical
    gsl_histogram_increment(histogram, value);
}

// incremene histogram with a soecific weight
void HistObservable::observe(const double& value, const double& weight)
{
    // increment histogram with weight
    #pragma omp critical
    gsl_histogram_accumulate(histogram, value, weight);
}

// returnr the value of bin i
double HistObservable::get_value(const int& i)
{
    // value to return
    double rval = 0.0;

    // get the value of bin i if less than number of bins
    if(i<N)
        rval = gsl_histogram_get(histogram, i);
    else
    {
        printf("ERROR: trying to access bin outside of range in Histogram!\n");
        exit(-1);
    }

    // retunr the value
    return rval;
}

// return the number of bins
int HistObservable::get_nbins(){return N;}

// return the center point of bin i
double HistObservable::get_bin_center(const int& i)
{
    double lower;
    double upper;
    double center = 0.0;

    if(i<N)
        gsl_histogram_get_range(histogram, i, &lower, &upper);
    else
    {
        printf("ERROR: trying to access bin center outside of range in Histogram!\n");
        exit(-1);
    }

    // calculate the center
    center = 0.5 * (upper + lower);

    return center;
}

double HistObservable::get_bin_width(const int& i)
{
    double lower;
    double upper;
    double dbin = 0.0;

    if(i<N)
        gsl_histogram_get_range(histogram, i, &lower, &upper);
    else
    {
        printf("ERROR: trying to access bin center outside of range in Histogram!\n");
        exit(-1);
    }

    // calculate the center
    dbin = (upper - lower);

    return dbin;
}

double HistObservable::get_pdf(const int& i)
{
    // get the bin width
    double dx = get_bin_width(i);

    // get the sum of the histogram
    double tsum = gsl_histogram_sum(histogram);

    // calculate the normalised bin height
    double pdf = this->get_value(i) / (tsum * dx);

    // return the normalised height
    return pdf;
}
