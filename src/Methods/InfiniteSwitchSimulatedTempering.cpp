#include "InfiniteSwitchSimulatedTempering.hpp"

using namespace::std;

/******************************************************************************
                    Infinite Switch Simulated Tempering
 *****************************************************************************/
InfiniteSwitchSimulatedTempering::InfiniteSwitchSimulatedTempering(
                                     Molecule* molecule_pt, const double& tmin,
                                     const double& tmax, const unsigned& nint,
                                     const double& t_step, const double& tau)
{
    // set the integration limits
    beta_max = 1.0 / tmin;
    beta_min = 1.0 / tmax;

    // initialise the Thermal Force Scaling
    thermal_force_scaling = 0.0;

    // set the number of integration points
    number_of_interpolation_points = nint;

    // set the alpha parameter
    time_step_scaling = tau;

    // set the timestep
    time_step = time_step_scaling * t_step;



    // get the legendre polynominals basis for the interpolation points
    gauss_weight = new double[number_of_interpolation_points];
    beta = new double[number_of_interpolation_points];
    partition_estimate = new AverageObservable[number_of_interpolation_points];

    legendre_compute_glr(number_of_interpolation_points, beta, gauss_weight);
    rescale(beta_min, beta_max, number_of_interpolation_points,
            beta, gauss_weight);

    // set the initial weights
    beta_weight = new double[number_of_interpolation_points];

    for(unsigned i=0; i<number_of_interpolation_points; i++)
        beta_weight[i] = 1.0;

    // rescale the temperature weights
    normalise_weights();

    // calculate the force rescaling
    calculate_force_rescaling(molecule_pt);

    // initialise the force
    init_force(molecule_pt);

    // update temperature
    molecule_pt->set_beta(beta[number_of_interpolation_points-1]);
}

// destructor
InfiniteSwitchSimulatedTempering::~InfiniteSwitchSimulatedTempering()
{
    for(unsigned i=0; i<number_of_interpolation_points; i++)
    {
        delete &beta[i];
        delete &gauss_weight[i];
        delete &beta_weight[i];
        delete &partition_estimate[i];
    }

    // double holders
    delete &beta_min;
    delete &beta_max;
    delete &time_step;
    delete &thermal_force_scaling;
    delete &time_step_scaling;

    // unsigned holders
    delete &number_of_interpolation_points;
}

// apply the isst force rescaling
void InfiniteSwitchSimulatedTempering::apply_force_rescaling(
                                                        Molecule* molecule_pt)
{
	// update the force rescaling using the molecule
    update_force_rescaling(molecule_pt);

	// dereference
	unsigned number_of_particles = molecule_pt->nparticle();
	unsigned dimension = molecule_pt->dim();
    Particle* particle_pt = NULL;

	#pragma omp simd
	for(unsigned j=0; j<number_of_particles; j++)
    {
        particle_pt = molecule_pt->particle_pt(j);
        // rescale the forces
        for(unsigned i=0; i<dimension; i++)
            *particle_pt->f_pt(i) *= thermal_force_scaling;

        // rescale the torque
        if(particle_pt->rigid_body() == true)
            *particle_pt->tau_pt(0) *= thermal_force_scaling;
    }

}

// get the i^th beta
double InfiniteSwitchSimulatedTempering::get_beta(const unsigned& i)
{
	return beta[i];
}

// get the observable weight used for reweighting observables
vector<double> InfiniteSwitchSimulatedTempering::get_observable_weights(
                                                        Molecule* molecule_pt)
{
	// parameters
  	double V = *molecule_pt->potential_pt();
	vector<double> Weights(number_of_interpolation_points, 0.0);

	// loop over integration point to find
   	for(unsigned i=0; i<number_of_interpolation_points; i++)
    {
        // dereference
		double z_i = partition_estimate[i].get_average();
       	double beta_i = beta[i];
		double integral = 0.0;

		// loop over integrals
		for(unsigned j=0;j<number_of_interpolation_points;j++)
		{
            // dereference
		    double w_j = beta_weight[j];
		    double g_j = gauss_weight[j];
		    double beta_j = beta[j];
		    // calculate integral
		    integral += g_j * w_j * exp(-(beta_j-beta_i) * V);
		}

       	// assign the weight
       	Weights[i] = 1.0 / (z_i * integral);
    }

    // return the Weight
    return Weights;
}

void InfiniteSwitchSimulatedTempering::set_beta_weight(
                                                const vector<double>& weight,
						                        Molecule* molecule_pt)
{
	// assign the weights
	for(unsigned i=0; i<number_of_interpolation_points; i++)
		beta_weight[i] = weight[i];

	// normalise the weights
	normalise_weights();
	// calculate the force rescaling
	calculate_force_rescaling(molecule_pt);
	// initialise the force
	init_force(molecule_pt);
}

// initialise the force rescaling
void InfiniteSwitchSimulatedTempering::calculate_force_rescaling(
                                                        Molecule* molecule_pt)
{
	// dereference
	double V = *molecule_pt->potential_pt();
	double InvTemp = molecule_pt->beta();

	double BarNumSum = 0.0;
	double BarDeNumSum = 0.0;

	for(unsigned i=0;i<number_of_interpolation_points;i++)
	{
	    // dereference the weight
	 	double w_i = beta_weight[i];

	    // sum up the integrals
	    BarNumSum   += gauss_weight[i] * beta[i] * w_i * exp(-beta[i] * V);
	    BarDeNumSum += gauss_weight[i] * w_i * exp(-beta[i] * V);
	 }

	 // update the thermal scaling
	 thermal_force_scaling = BarNumSum / (InvTemp * BarDeNumSum);
};

// initialise the force
void InfiniteSwitchSimulatedTempering::init_force(Molecule* molecule_pt)
{
  unsigned nParticle = molecule_pt->nparticle();
  unsigned Dim = molecule_pt->dim();

  for(unsigned j=0;j<nParticle;j++)
    for(unsigned i=0;i<Dim;i++)
      *molecule_pt->particle_pt(j)->f_pt(i) *= thermal_force_scaling;
};

// learn the beta weights
void InfiniteSwitchSimulatedTempering::learn_beta_weights(Molecule* molecule_pt)
{
	// dereference
	double h = time_step;
	double V = *molecule_pt->potential_pt();
	double BarDeNumSum = 0.0;

	for(unsigned i=0; i<number_of_interpolation_points; i++)
	{
		// dereference the weight
		double w_i = beta_weight[i];

		// sum up the integrals
		BarDeNumSum += gauss_weight[i] * w_i * exp(-beta[i] * V);
	}

	// update the partition estimate
	for(unsigned i=0; i<number_of_interpolation_points; i++)
	{
		// partition estimate
		partition_estimate[i].observe(exp(-beta[i] * V) / BarDeNumSum);
	}

	// update all the weights with the new values
	for(unsigned i=0; i<number_of_interpolation_points; i++)
	{
		// new weights
		double omega_np1 = 0.0;

		// current weight
		double omega_n = beta_weight[i];

		// get the partion partition
		double z_n = partition_estimate[i].get_average();

		// calculate the new f's
		omega_np1 = (1.0 - h) * omega_n + h / z_n;

		// update the new weights
		beta_weight[i] = omega_np1;
	}

	// normalise the weights
	normalise_weights();
}

// normalise the beta weights
void InfiniteSwitchSimulatedTempering::normalise_weights()
{
	// weight normalisation
	double weight_norm = 0.0;

	// integrate over the weights
	for(unsigned i=0; i<number_of_interpolation_points; i++)
		weight_norm += gauss_weight[i] * beta_weight[i];

	// normalise all the weights
	for(unsigned i=0; i<number_of_interpolation_points; i++)
		beta_weight[i] /= weight_norm;
}

// update the force rescaling
void InfiniteSwitchSimulatedTempering::update_force_rescaling(
                                                    Molecule* molecule_pt)
{
	// learn the beta weights
	learn_beta_weights(molecule_pt);

	// calculate the force rescaling
	calculate_force_rescaling(molecule_pt);
};
