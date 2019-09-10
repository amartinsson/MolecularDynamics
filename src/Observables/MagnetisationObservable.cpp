#include "MagnetisationObservable.hpp"

// constructor
MagnetisationObservable::MagnetisationObservable(Molecule* molecule_pt,
    const double& mmin, const double& mmax, const int& N, const int& recf,
        const int& rect)
            : SystemObservable(recf, rect),
                magDist(HistObservable(mmin, mmax, N)), molecule_pt(molecule_pt)
{};
// update function
void MagnetisationObservable::update()
{
    if(recStep()) {
        magDist.observe(molecule_pt->magnetisation());
    }
}

// update function
void MagnetisationObservable::update(const double& weight)
{
    if(recStep()) {
        magDist.observe(molecule_pt->magnetisation(), weight);
    }
}

// print function
void MagnetisationObservable::print(const char* file_name,
    const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    // get the number of bins
    int N = magDist.get_nbins();

    for(int i=0; i<N; i++) {
        fprintf(file, "%1.7e, %1.7e, %1.7e\n",
                magDist.get_bin_center(i),
                magDist.get_value(i),
                magDist.get_pdf(i));
    }

    // close the file
    fclose(file);
}

// get average
double MagnetisationObservable::get_average()
{
    printf("ERROR: magnetisation cannot use get_average!\n");
    exit(-1);
}
// get instant
double MagnetisationObservable::get_instant()
{
    printf("ERROR: magnetisation cannot use get_instant!\n");
    exit(-1);
}
