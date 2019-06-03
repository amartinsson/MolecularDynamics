#include "SystemTemperature.hpp"

// constructor
SystemTemperature::SystemTemperature(Molecule* molecule_pt) :
    momentum(false), configurational(false)
{
    // copy the memory address of the molecule
    SystemTemperature::system = molecule_pt;

    // momentum temperature
    mTemp = new AverageObservable();

    // configurational temperature
    cTemp = new AverageObservable();

}

// get instant function
double SystemTemperature::get_instant()
{
    double instant = 0.0;

    if(momentum)
        instant = mTemp->get_instant();
    else if(configurational)
        instant = cTemp->get_instant();
    else
    {
        printf("ERROR: Need to specify which temperature to calulate!\n");
        exit(-1);
    }

    return instant;
}

// get instant function
double SystemTemperature::get_average()
{
    double average = 0.0;

    if(momentum)
        average = mTemp->get_average();
    else if(configurational)
        average = -cTemp->get_average();
    else
    {
        printf("ERROR: Need to specify which temperature to calulate!\n");
        exit(-1);
    }

    return average;
}


// update function
void SystemTemperature::update()
{
    if(momentum)
        update_mtemp();
    else if(configurational)
        update_ctemp();
    else
    {
        printf("ERROR: Need to specify which temperature to calulate!\n");
        exit(-1);
    }
}

// print function at timeindex
void SystemTemperature::print(const char* file_name, const double& time,
                              const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* temperature_file = fopen(filename, "a");

    fprintf(temperature_file, "%1.7e, %1.7e, %1.7e\n",
            time, this->get_instant(), this->get_average());

    // close the file
    fclose(temperature_file);
}

// print function at time
void SystemTemperature::print(const char* file_name, const double& time)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s.csv", (name).c_str());
    FILE* temperature_file = fopen(filename, "a");

    fprintf(temperature_file, "%1.7e, %1.7e, %1.7e\n",
            time, this->get_instant(), this->get_average());

    // close the file
    fclose(temperature_file);

}

// calculate momentum temperature
void SystemTemperature::update_mtemp()
{
    // make an average observable
    AverageObservable temp = AverageObservable();

    // loop overa all the particles
    for(const auto& particle : system->Particles)
        temp.observe(particle.second->p.dot(
                     particle.second->m.inv() * particle.second->p));

    // add the observed temperature to the mTemp observable
    double dim = double(system->dim());
    double N = double(system->nparticle());

    mTemp->observe(N * temp.get_average() / (dim * (N - 1.0)));
}

// calculate configurational temperature
void SystemTemperature::update_ctemp()
{
    // make an average observable
    AverageObservable temp = AverageObservable();

    // loop overa all the particles
    for(const auto& particle : system->Particles)
        temp.observe(particle.second->q.dot(particle.second->f));

    // add the observed temperature to the cTemp observable
    // double dim = double(system->dim());
    // cTemp->observe(temp.get_average() / dim);
    cTemp->observe(temp.get_average());
}
