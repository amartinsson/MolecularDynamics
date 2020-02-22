#include "SphereOrderObservable.hpp"

using namespace::std;

// constructor
SphereOrderObservable::SphereOrderObservable(Molecule* molecule_pt,
                                             const int& l_max,
                                             const double& coff,
                                             const int& recf, const int& rect)
    : SystemObservable(recf, rect), lmax(l_max), cutoff(coff),
        looplmax(2 * l_max + 1), locRecStep(false)
{
    // make array of particles
    SphereOrderObservable::system = molecule_pt;

    // allocate memory in unorder maps structure
    this->allocate_unordered_maps();

    // make a new average state over all particles
    state = new AverageObservable();
}

SphereOrderObservable::~SphereOrderObservable()
{
    // remove data structure
    for(const auto& particle : system->Particles)
    {
        for(auto& ave : obsReal.at(particle.second))
            delete ave;

        for(auto& ave : obsImag.at(particle.second))
            delete ave;

        obsReal.at(particle.second).clear();
        obsReal.at(particle.second).clear();
    }

    delete state;
}

// clear the observations
void SphereOrderObservable::clear()
{
    // clear all the observales
    for(const auto& particle : system->Particles)
        for(unsigned i=0; i<looplmax; i++)
        {
            obsReal.at(particle.second)[i]->clear();
            obsImag.at(particle.second)[i]->clear();
        }

    // call recStep to track if to update
    locRecStep = recStep();
}

// update function
void SphereOrderObservable::update()
{
    printf("ERROR: must update with value for r in SphereOrderObservable\n");
    exit(-1);
}

void SphereOrderObservable::update(Particle* current, Particle* neighbour,
                                   const double& r, const Vector& dr)
{
    if(r < cutoff && locRecStep)
    {
        // calculate the angles between particles
        vector<double> angles = get_angles(dr.unit());

        // add observable to current particle
        add_obs(current, angles);

        // reverse direction for neighbour particle
        angles = get_angles((dr.neg()).unit());

        // add observable to neighbour
        add_obs(current, angles);
    }
}

// print function at time_index
void SphereOrderObservable::print(const char* file_name, const double& time,
                                  const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // // calculate average for all the particles
    AverageObservable current_state = AverageObservable();

    // loop over all the particles
    for(const auto& particle : system->Particles)
    {
        // get the real an imaginary part
        double real = 0.0;
        double imag = 0.0;

        for(unsigned i=0; i<looplmax; i++)
        {
            real += obsReal.at(particle.second)[i]->get_average();
            imag += obsImag.at(particle.second)[i]->get_average();
        }

        // calculate the absolute value
        double obs = std::sqrt(4.0 * M_PI / (2.0 * (double)lmax + 1.0) * (real * real + imag * imag));

        // reord the observation
        current_state.observe(obs);
    }

    // add to the state average observable
    state->observe(current_state.get_average());

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    // print the current state and long term average
    fprintf(file, "%1.7e, %1.7e, %1.7e\n", time, state->get_instant(),
            state->get_average());

    // close the file
    fclose(file);
}

// print function at time
void SphereOrderObservable::print(const char* file_name, const double& time)
{
    printf("ERROR: use other print funtion!\n");
    exit(-1);
}

// get the instant order value
double SphereOrderObservable::get_instant()
{
    return state->get_instant();
}

// get the average order value
double SphereOrderObservable::get_average()
{
    return state->get_average();
}

// allocate memory for Real and Imaginary parts
void SphereOrderObservable::allocate_unordered_maps()
{
    // reserve memory in unordered map
    obsReal.reserve(system->nparticle());
    obsImag.reserve(system->nparticle());

    // add particles to memory
    for(const auto& particle : system->Particles)
    {
        // make vector of length we need
        vector<AverageObservable*> vec(looplmax);

        // assign new memory
        for(auto& el : vec)
            el = new AverageObservable();

        // insert into the unordered map
        obsReal.insert(make_pair(particle.second, vec));

        // assign new memory
        for(auto& el : vec)
            el = new AverageObservable();

        // insert into the unordered map
        obsImag.insert(make_pair(particle.second, vec));
    }
}

// physics is stupid, switched order of theta and phi!!! check wikipedia
vector<double> SphereOrderObservable::get_angles(const Vector& dr)
{
    // store for angles
    vector<double> angles(2, 0.0);

    // calculate theta first
    if(dr.size() > 2)
        angles[0] = acos(dr(2));
    else
        angles[0] = 0.0;

    // calculate phi second
    angles[1] = atan2(dr(1), dr(0));

    return angles;
}

// add an observation to the particle with angles
void SphereOrderObservable::add_obs(Particle* particle,
                                    const vector<double>& angles)
{
    // calculate the order parameter
    Matrix obs = calculate_order_param(angles);

    for(unsigned i=0; i<looplmax; i++)
    {
        obsReal.at(particle)[i]->observe(obs(i, 0));
        obsImag.at(particle)[i]->observe(obs(i, 1));
    }
}

Matrix SphereOrderObservable::calculate_order_param(const vector<double>& angles)
{
    // dereference some helpers
    double theta = angles[0];
    double phi = angles[1];

    // make matrix to return
    Matrix RetMat(looplmax, 2);

    unsigned index = 0;

    for(unsigned m=0; m<=lmax; m++)
    {
        // make the spherical harmonic
        complex<double> sph = boost::math::spherical_harmonic(lmax, m,
                                                              theta, phi);

        RetMat(index, 0) = real(sph);
        RetMat(index, 1) = imag(sph);
         if(m > 0)
         {
             // bump index
             index++;

             // calculate real and imaginary parts
             RetMat(index, 0) = pow(-1.0, m) * real(sph);
             RetMat(index, 1) = pow(-1.0, m + 1) * imag(sph);
        }

        // bump index
        index++;
    }

    return RetMat;
}

// print function at time_index
void SphereOrderObservable::print_traj(const char* file_name, const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/Frames/%s_%i.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "w");

    // loop over all the particles
    for(const auto& particle : system->Particles)
    {
        // get the real an imaginary part
        double real = 0.0;
        double imag = 0.0;

        for(unsigned i=0; i<looplmax; i++)
        {
            real += obsReal.at(particle.second)[i]->get_average();
            imag += obsImag.at(particle.second)[i]->get_average();
        }

        // calculate the absolute value
        double obs = std::sqrt(4.0 * M_PI / (2.0 * (double)lmax + 1.0) * (real * real + imag * imag));

        // record the observation
        fprintf(file, "%.4f\n", obs);
    }

    // close the file
    fclose(file);
}
