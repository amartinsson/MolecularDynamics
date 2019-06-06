#include "SystemObservable.hpp"

using namespace::std;

void SystemObservable::print(const char* file_name, const double& time,
                             const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    fprintf(file, "%1.7e, %1.7e, %1.7e\n",
            time, get_instant(), get_average());

    // close the file
    fclose(file);
}

void SystemObservable::print(const char* file_name, const double& time)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s.csv", (name).c_str());
    FILE* file = fopen(filename, "a");

    fprintf(file, "%1.7e, %1.7e, %1.7e\n",
            time, this->get_instant(), this->get_average());

    // close the file
    fclose(file);
}

// function to check that we need to record this step
bool SystemObservable::recStep()
{
    // boolean
    bool Record = false;

    // check if this is a record step
    if(RecCount % RecFreq == 0 && RecCount != 0)
    {
        // set record to true
        Record = true;

        // reset the counter
        RecCount = 1;

        // check that we are above threshold
        if(RecTotCount < RecThresh)
            Record = false;
    }
    else
        RecCount++;

    // increment the total counter
    if(RecTotCount < RecThresh)
        RecTotCount++;

    // return the correct value
    return Record;
}
