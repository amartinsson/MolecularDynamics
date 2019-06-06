#include "OrderObservable.hpp"

using namespace::std;

// constructor
OrderObservable::OrderObservable(Molecule* molecule_pt, const double& part_rad,
                                 const int& recf, const int& rect)
    : SystemObservable(recf, rect), cutoff((1.0 + std::sqrt(2.0)) * part_rad),
        pr(0.5 * part_rad)
{
    // copy the pointer to the moleucle object
    OrderObservable::system = molecule_pt;

    // reserve the number of positions in memory
    observablemapReal.reserve(system->nparticle());
    observablemapImag.reserve(system->nparticle());

    // make map of particles
    for(const auto& particle : system->Particles)
    {
        observablemapReal.insert(make_pair(particle.second,
                                       new AverageObservable()));
        observablemapImag.insert(make_pair(particle.second,
                                       new AverageObservable()));
    }

    // set wavevector
    set_wave();
    // initialise state vector
    state = new AverageObservable();
}

// destructor
OrderObservable::~OrderObservable()
{
    for(const auto& particle : system->Particles)
    {
        delete observablemapReal.at(particle.second);
        delete observablemapImag.at(particle.second);
    }

    delete state;
}

// clear the observations
void OrderObservable::clear()
{
    for(const auto& particle : system->Particles)
    {
        observablemapReal.at(particle.second)->clear();
        observablemapImag.at(particle.second)->clear();
    }

    // bump recstep
    localRecStep = recStep();
}

// wrong update function
void OrderObservable::update()
{
    printf("ERROR: must update with value for r OrderObservable\n");
    exit(-1);
}

void OrderObservable::update(Particle* current, Particle* neighbour,
                             const double& r, const Vector& dr)
{
    // check that we are smaller than the cutoff
    if(r < cutoff && localRecStep)
    {
        // get the order parameter for separation
        vector<double> order = get_orderparam(dr);

        // add observation to current particle
        observablemapReal.at(current)->observe(order[0]);
        observablemapImag.at(current)->observe(order[1]);

        // add observation to neighbour particle
        observablemapReal.at(neighbour)->observe(order[0]);
        observablemapImag.at(neighbour)->observe(-order[1]);
    }
}

// get the observable for particular distance
vector<double> OrderObservable::get_orderparam(const Vector& dr)
{
    // return value
    vector<double> retval(2, 0.0);

    // loop over all the vectors
    for(unsigned i=0; i<6; i++)
    {
        double theta = k[i].dot(dr);

        retval[0] += 1.0 / 6.0 * cos(theta);
        retval[1] += 1.0 / 6.0 * sin(theta);
    }

    // return
    return retval;
}


// print function at time_index
void OrderObservable::print(const char* file_name, const double& time,
                            const unsigned& index)
{
    // convert the filename to a string
    std::string name(file_name);

    // open the file to write to
    char filename[50];
    sprintf(filename, "Observables/States/%s_%d.csv", (name).c_str(), index);
    FILE* file = fopen(filename, "a");

    // calculate average for all the particles
    AverageObservable current_state = AverageObservable();

    // loop over all the particles
    for(const auto& particle : system->Particles)
    {
        // get the real an imaginary part
        double real = observablemapReal.at(particle.second)->get_average();
        double imag = observablemapImag.at(particle.second)->get_average();

        // calculate the absolute value
        double obs = real * real + imag * imag;

        // reord the observation
        current_state.observe(obs);

        fprintf(file, "%1.7e\n", obs);
    }

    fclose(file);

    // add to the state average observable
    state->observe(current_state.get_average());

    // open the file to write to
    // char filename[50];
    sprintf(filename, "Observables/%s_0.csv", (name).c_str());
    // FILE* file = fopen(filename, "a");
    file = fopen(filename, "a");

    // print the current state and long term average
    fprintf(file, "%1.7e, %1.7e, %1.7e\n", time, state->get_instant(),
            state->get_average());

    // close the file
    fclose(file);
}

// print function at time
void OrderObservable::print(const char* file_name, const double& time)
{
    printf("ERROR: use other print funtion!\n");
    exit(-1);
}

// get the instant order value
double OrderObservable::get_instant()
{
    return state->get_instant();
}

// get the average order value
double OrderObservable::get_average()
{
    return state->get_average();
}

// set the wave numbers
void OrderObservable::set_wave()
{
    // resize and initialise the vector
    k.resize(6, Vector(3));

    // wave vector 0
    k[0](0) = -M_PI / pr;
    k[0](1) = M_PI / pr;
    k[0](2) = 0.0;
    // wave vector 1
    k[1](0) = -M_PI / pr;
    k[1](1) = 0.0;
    k[1](2) = M_PI / pr;
    // wave vector 2
    k[2](0) = 0.0;
    k[2](1) = M_PI / pr;
    k[2](2) = M_PI / pr;
    // wave vector 3
    k[3](0) = -(3.0 * M_PI / (2.0 * pr));
    k[3](1) = 3.0 * M_PI / (2.0 * pr);
    k[3](2) = M_PI / (2.0 * pr);
    // wave vector 4
    k[4](0) = -(3.0 * M_PI / (2.0 * pr));
    k[4](1) = M_PI / (2.0 * pr);
    k[4](2) = 3.0 * M_PI / (2.0 * pr);
    // wave vector 5
    k[5](0) = -(M_PI / (2.0 * pr));
    k[5](1) = 3.0 * M_PI / (2.0 * pr);
    k[5](2) = 3.0 * M_PI / (2.0 * pr);
}
