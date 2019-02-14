#include "SimulatedTempering.hpp"

using namespace::std;

/******************************************************************************
                            Simulated Tempering
 *****************************************************************************/
SimulatedTempering::SimulatedTempering(const double& tmin, const double& tmax,
				   					   const double& n_temperatures,
				   		   			   const unsigned& mod_switch,
								   	   const int& seed) :
				                       uniform_gen(0.0, 1.0, seed)
{
	// read in all the parameteres
	beta_min = 1.0 / tmax;
	beta_max = 1.0 / tmin;

	number_of_temperatures = n_temperatures;
	switch_frequency = mod_switch;

	// set current temperature index to the lowest temperature
	current_index = number_of_temperatures - 1;

	// set all the temperatures
	double delta_beta = (beta_max - beta_min) / ((double)n_temperatures - 1.0);

	// assign all the temperatures
	for(unsigned i=0; i<number_of_temperatures; i++)
	{
		// assigne weights and temperatures
		beta.push_back(beta_min + (double)i * delta_beta);
		Weights.push_back(1.0);

		// assigne the number of steps counter
		state_visited.push_back(false);
		ave_V.push_back(new AverageObservable);
	}

	// initialise random number generator
	// long unsigned seed;
	// struct timeval tv;
	// gettimeofday(&tv,0);
	// seed = tv.tv_sec + tv.tv_usec + 7 * omp_get_thread_num();
	//
	// uniform = gsl_rng_alloc(gsl_rng_mt19937);
	// gsl_rng_set(uniform, seed);
}

// destructor
SimulatedTempering::~SimulatedTempering()
{
    // free memory for vectors
    for(unsigned i=0; i< number_of_temperatures; i++)
    {
        delete &beta[i];
        delete &Weights[i];
        delete ave_V[i];
    }
    beta.clear();
    Weights.clear();
    state_visited.clear();
    ave_V.clear();

    // free memory for doubles and ints
    delete &beta_min;
    delete &beta_max;
    delete &current_index;
    delete &number_of_temperatures;
    delete &switch_frequency;

    // free uniform random
	delete &uniform_gen;
    // gsl_rng_free(uniform);
}

// get the i^th beta
double SimulatedTempering::get_beta(const unsigned& i)
{
	return beta[i];
}

// return the temperature index of the current temperature
unsigned SimulatedTempering::get_temperature_index()
{
	return current_index;
}

// return the current temperature of the simulation
double SimulatedTempering::get_temperature()
{
	return 1.0 / beta[current_index];
}

// set the temperature weight for i^th component
void SimulatedTempering::set_beta_weight(const double& weight,
										 const unsigned& i)
{
	Weights[i] = weight;
}

// normalise the weights
void SimulatedTempering::normalise_weights()
{
	// dereference the constant
    double weight_sum = 0.0;

	// sum over the temperatures
    for(unsigned i=0; i<number_of_temperatures; i++)
		weight_sum += Weights[i];

    // set the rescaled weights
    for(unsigned i=0; i<number_of_temperatures; i++)
		Weights[i] /= weight_sum;
}

void SimulatedTempering::update_temperature(Molecule* molecule_pt,
											const unsigned& step)
{
	// if we're on a switch step then swith
	// otherwise do nothing
	if(step % switch_frequency == 0)
	{
		// make proposed index
		double V = *molecule_pt->potential_pt();
		double U = 0.0;
		int shift = 0;

		// check if the current step has been visited
		if(state_visited[current_index] == false)
			state_visited[current_index] = true;

		// observe the potential energy for this state
		ave_V[current_index]->observe(V);

		// calculate the weights for this step
		for(unsigned i=1; i<number_of_temperatures; i++)
		{
			// If a state has not been visited, assume it has the same PE as preceding state
            // to increase transition probability (see http://dx.doi.org/10.1021/acs.jpcb.5b03381 )
			if(state_visited[i] != false)
				Weights[i] = Weights[i-1] * exp(0.5 * (beta[i] - beta[i-1])
					   * (ave_V[i]->get_average() + ave_V[i-1]->get_average()));
			else
				Weights[i] = Weights[i-1] * exp(0.5 * (beta[i] - beta[i-1])
												* ave_V[i-1]->get_average());
		}

		// generator uniform random variable
		U = uniform_gen();//gsl_rng_uniform(uniform);

		// propose shift up
		if(U >= 0.5 && current_index + 1 < number_of_temperatures)
			shift = 1;

		// propose shift down
		if(U < 0.5 && current_index > 0)
			shift = -1;

		// if the shift is non-zero
		if(shift != 0)
		{
			// calculate acceptance probability
			double accept = Weights[current_index + shift]
							/ Weights[current_index]
				* exp(-(beta[current_index + shift] - beta[current_index]) * V);

			 // if we accept the state, update the temperatures
		  	// and the index - else do nothing
		  	if(U <= accept)
		  	{
		  		// update temperature index
		  		current_index += shift;

				// rescale the momentum
				rescale_momentum(molecule_pt, shift);

		  		// update temperature
		  		molecule_pt->set_beta(beta[current_index]);
		  	}
		}
	}
}

void SimulatedTempering::rescale_momentum(Molecule* molecule_pt, const int& shift)
{
	// get the previous index
	int previous_index = current_index - shift;
	unsigned number_of_particles = molecule_pt->nparticle();
	unsigned dim = molecule_pt->dim();
	double* p = NULL;

	// calculate the scale
	double scale = sqrt(beta[previous_index] / beta[current_index]);

	Particle* particle_i = NULL;

	for(unsigned i=0; i<number_of_particles; i++)
		for(unsigned j=0; j<dim; j++)
		{
			particle_i = molecule_pt->particle_pt(i);
			p = particle_i->p_pt(j);

			(*p) *= scale;
		}
}
