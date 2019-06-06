#include "RadialDistObservable.hpp"

using namespace::std;

RadialDistObservable::RadialDistObservable(const double& rmin,
                                           const double& rmax,
                                           const int& N, const int& recf,
                                           const int& rect)
                                : SystemObservable(recf, rect),
                                    radialDist(HistObservable(rmin, rmax, N))
{
    // initialise the histogram observable with the correct parameters
}

void RadialDistObservable::bump_recstep()
{
    localRecStep = recStep();
}


void RadialDistObservable::update()
{
    printf("ERROR: must update with value for r in radial RadialDistObservable\n");
    exit(-1);
}


void RadialDistObservable::update(const double& r)
{
    if(localRecStep)
        // observe the separation
        radialDist.observe(r);
}

void RadialDistObservable::print(const char* file_name, const double& time,
                                 const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    // get the number of bins
    int N = radialDist.get_nbins();

    for(int i=0; i<N; i++)
        fprintf(file, "%1.7e, %1.7e %1.7e\n",
                radialDist.get_bin_center(i),
                radialDist.get_value(i),
                radialDist.get_pdf(i));

    // close the file
    fclose(file);
}

// get average
double RadialDistObservable::get_average()
{
    printf("ERROR: radial cannot use get_average!\n");
    exit(-1);
}
// get instant
double RadialDistObservable::get_instant()
{
    printf("ERROR: radial cannot use get_instant!\n");
    exit(-1);
}
