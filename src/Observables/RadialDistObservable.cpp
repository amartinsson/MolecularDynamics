#include "RadialDistObservable.hpp"

using namespace::std;

RadialDistObservable::RadialDistObservable(const double& rmin,
                                           const double& rmax,
                                           const int& N, const int& recf,
                                           const int& rect,
                                           const int& nparticles)
                                : SystemObservable(recf, rect),
                                    radialDist(HistObservable(rmin, rmax, N)),
                                        number_of_particles(nparticles)
{
    // initialise the average number density holder
    number_density = new AverageObservable();
}

void RadialDistObservable::bump_recstep(const double& V)
{
    // bump the recStep
    localRecStep = recStep();

    // observe the number density
    if(localRecStep)
        number_density->observe((double)number_of_particles / V);
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
        fprintf(file, "%1.7e, %1.7e, %1.7e\n",
                radialDist.get_bin_center(i),
                radialDist.get_value(i),
                this->get_rpdf_3d(i));

    // close the file
    fclose(file);
}

double RadialDistObservable::get_rpdf_3d(const int& i)
{
    // get the bin count
    double bin_count = radialDist.get_value(i);
    // get average number density
    double rho = number_density->get_average();
    // double rho = number_density->get_instant();
    // get the lower bin
    double r = radialDist.get_lower_bin(i);
    // get the bn width
    double dr = radialDist.get_bin_width(i);

    // return the radial distribution normalised
    return radialDist.get_pdf(i) / (rho * 4.0 * M_PI * pow(r, 2.0) * dr);
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
