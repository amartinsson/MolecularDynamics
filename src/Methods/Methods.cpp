// // Include the methods Header file
//
// #include "Methods.hpp"
// #include <algorithm>
// 
// // ---------------------------------------------------------------------------//
// //                     Infinite Swap Simulated Tempering
// // ---------------------------------------------------------------------------//
//
// // constructor for the ISST method
// ISST::ISST(Molecule* MolPt, const double& tmin, const double& tmax,
// 	       const double& nInt, const double& t_step, const double& Tau)
// {
//   // set the integration limits
//   beta_max = 1.0 / tmin;
//   beta_min = 1.0 / tmax;
//
//   // initialise the Thermal Force Scaling
//   thermal_force_scaling = 0.0;
//
//   // set the number of integration points
//   number_of_interpolation_points = nInt;
//
//   // set the alpha parameter
//   time_step_scaling = Tau;
//
//   // set the timestep
//   time_step = time_step_scaling * t_step;
//
//   // make the integration table
//   gsl_integration_glfixed_table * tbl =
//   			gsl_integration_glfixed_table_alloc(number_of_interpolation_points);
//
//   // set the integration points and the reciprocal temperatures
//   for(unsigned i=0; i<number_of_interpolation_points; i++)
//     {
//       // helpers
//       double x_i = 0.0;
//       double w_i = 0.0;
//
//       // get the integration points
//       gsl_integration_glfixed_point(beta_min, beta_max, i, &x_i, &w_i, tbl);
//
//       // add the new values to vector
//       beta.push_back(x_i);
//       gauss_weight.push_back(w_i);
//
//       // set the temperature weights
// 	  beta_weight.push_back(1.0);
//
// 	  // initialise the average partion Observables
// 	  partition_estimate.push_back(new AverageObservable());
//     }
//
//   // free the integration table
//   gsl_integration_glfixed_table_free(tbl);
//
//   // rescale the temperature weights
//   normalise_weights();
//
//   // calculate the force rescaling
//   calculate_force_rescaling(MolPt);
//
//   // initialise the force
//   init_force(MolPt);
//
//   // update temperature
//   MolPt->set_beta(beta[number_of_interpolation_points-1]);
//   //MolPt->set_beta(beta[0]);
// };
//
// // apply the isst force rescaling
// void ISST::apply_force_rescaling(Molecule* molecule_pt)
// {
// 	// update the force rescaling using the molecule
//     update_force_rescaling(molecule_pt);
//
// 	// dereference
// 	unsigned number_of_particles = molecule_pt->nparticle();
// 	unsigned dimension = molecule_pt->dim();
//
// 	// rescale the force
// 	#pragma omp simd collapse(2)
// 	for(unsigned j=0; j<number_of_particles; j++)
//     	for(unsigned i=0; i<dimension; i++)
//   			*molecule_pt->particle_pt(j)->f_pt(i) *= thermal_force_scaling;
// }
//
// // dummy pure virtual function
// void ISST::get_method_info()
// {
// #pragma omp master
// {
// 	printf("-Accelerating integration with ISST:\n");
//
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 		printf("\tTemperature[%2d] = %1.2f, %1.7f\n", i, 1.0/beta[i], beta[i]);
//
// 	printf("--------------------------------------\n");
// }
// };
//
// void ISST::init_force(Molecule* MolPt)
// {
//   unsigned nParticle = MolPt->nparticle();
//   unsigned Dim = MolPt->dim();
//
//   for(unsigned j=0;j<nParticle;j++)
//     for(unsigned i=0;i<Dim;i++)
//       *MolPt->particle_pt(j)->f_pt(i) *= thermal_force_scaling;
// };
//
// // initialise the force rescaling
// void ISST::calculate_force_rescaling(Molecule* MolPt)
// {
// 	// dereference
// 	double V = *MolPt->potential_pt();
// 	double InvTemp = MolPt->beta();
//
// 	double BarNumSum = 0.0;
// 	double BarDeNumSum = 0.0;
//
// 	for(unsigned i=0;i<number_of_interpolation_points;i++)
// 	{
// 	    // dereference the weight
// 	 	double w_i = beta_weight[i];
//
// 	    // sum up the integrals
// 	    BarNumSum   += gauss_weight[i] * beta[i] * w_i * exp(-beta[i] * V);
// 	    BarDeNumSum += gauss_weight[i] * w_i * exp(-beta[i] * V);
// 	 }
//
// 	 // update the thermal scaling
// 	 thermal_force_scaling = BarNumSum / (InvTemp * BarDeNumSum);
// };
//
// // update the force rescaling
// void ISST::update_force_rescaling(Molecule* molecule_pt)
// {
// 	// learn the beta weights
// 	learn_beta_weights(molecule_pt);
//
// 	// calculate the force rescaling
// 	calculate_force_rescaling(molecule_pt);
// };
//
// void ISST::learn_beta_weights(Molecule* molecule_pt)
// {
// 	// dereference
// 	double h = time_step;
// 	double V = *molecule_pt->potential_pt();
// 	double BarDeNumSum = 0.0;
//
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 	{
// 		// dereference the weight
// 		double w_i = beta_weight[i];
//
// 		// sum up the integrals
// 		BarDeNumSum += gauss_weight[i] * w_i * exp(-beta[i] * V);
// 	}
//
// 	// update the partition estimate
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 	{
// 		// partition estimate
// 		partition_estimate[i]->observe(exp(-beta[i] * V) / BarDeNumSum);
// 	}
//
// 	// update all the weights with the new values
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 	{
// 		// new weights
// 		double omega_np1 = 0.0;
//
// 		// current weight
// 		double omega_n = beta_weight[i];
//
// 		// get the partion partition
// 		double z_n = partition_estimate[i]->get_average();
//
// 		// calculate the new f's
// 		omega_np1 = (1.0 - h) * omega_n + h / z_n;
//
// 		// update the new weights
// 		beta_weight[i] = omega_np1;
// 	}
//
// 	// normalise the weights
// 	normalise_weights();
// }
//
// std::vector<double> ISST::get_observable_weights(Molecule* MolPt)
// {
// 	// parameters
//   	double V = *MolPt->potential_pt(); //senare;
// 	std::vector<double> Weights(number_of_interpolation_points, 0.0);
//
// 	// loop over integration point to find
//    	for(unsigned i=0; i<number_of_interpolation_points; i++)
//     {
// 		// dereference
// 		double z_i = partition_estimate[i]->get_average();
//        	double beta_i = beta[i];
//
// 		double integral = 0.0;
//
// 		// loop over integrals
// 		for(unsigned j=0;j<number_of_interpolation_points;j++)
// 		{
// 		   // dereference
// 		   double w_j = beta_weight[j];
// 		   double g_j = gauss_weight[j];
// 		   double beta_j = beta[j];
//
// 		   // calculate integral
// 		   integral += g_j * w_j * exp(-(beta_j-beta_i) * V);
// 		}
//
//        	// assign the weight
//        	Weights[i] = 1.0 / (z_i * integral);
//      }
//
//    // return the Weight
//    return Weights;
// }
//
// // normalise the beta weights
// void ISST::normalise_weights()
// {
// 	// weight normalisation
// 	double weight_norm = 0.0;
//
// 	// integrate over the weights
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 		weight_norm += gauss_weight[i] * beta_weight[i];
//
// 	// normalise all the weights
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 		beta_weight[i] /= weight_norm;
// }
//
// // get the i^th beta
// double ISST::get_beta(const unsigned& i)
// {
// 	return beta[i];
// }
//
// // set the temperature weight for i^th component
// void ISST::set_beta_weight(const std::vector<double>& weight,
// 						   Molecule* molecule_pt)
// {
// 	// assign the weights
// 	for(unsigned i=0; i<number_of_interpolation_points; i++)
// 		beta_weight[i] = weight[i];
//
// 	// normalise the weights
// 	normalise_weights();
//
// 	// calculate the force rescaling
// 	calculate_force_rescaling(molecule_pt);
//
// 	// initialise the force
// 	init_force(molecule_pt);
// }
//
// // ---------------------------------------------------------------------------//
// //                              Simulated Tempering
// // ---------------------------------------------------------------------------//
//
// void SimulatedTempering::get_method_info()
// {
// 	#pragma omp master
// 	{
// 		printf("-Accelerating integration with Simulated Tempering:\n");
//
// 		for(unsigned i=0; i<number_of_temperatures; i++)
// 			printf("\tTemperature[%2d] = %1.2f\n", i, 1.0/beta[i]);
//
// 		printf("-----------------------------------------------------\n");
// 	}
// }
//
// // get the i^th beta
// double SimulatedTempering::get_beta(const unsigned& i)
// {
// 	return beta[i];
// }
//
// // return the temperature index of the current temperature
// unsigned SimulatedTempering::get_temperature_index()
// {
// 	return current_index;
// }
//
// // return the current temperature of the simulation
// double SimulatedTempering::get_temperature()
// {
// 	return 1.0 / beta[current_index];
// }
//
// // set the temperature weight for i^th component
// void SimulatedTempering::set_beta_weight(const double& weight,
// 										 const unsigned& i)
// {
// 	Weights[i] = weight;
// }
//
// // constructor
// SimulatedTempering::SimulatedTempering(const double& tmin, const double& tmax,
// 				   					   const double& n_temperatures,
// 				   		   			   const unsigned& mod_switch)
// {
// 	// read in all the parameteres
// 	beta_min = 1.0 / tmax;
// 	beta_max = 1.0 / tmin;
//
// 	number_of_temperatures = n_temperatures;
// 	switch_frequency = mod_switch;
//
// 	// set current temperature index to the lowest temperature
// 	current_index = number_of_temperatures - 1;
//
// 	// set all the temperatures
// 	double delta_beta = (beta_max - beta_min) / ((double)n_temperatures - 1.0);
//
// 	// assign all the temperatures
// 	for(unsigned i=0; i<number_of_temperatures; i++)
// 	{
// 		// assigne weights and temperatures
// 		beta.push_back(beta_min + (double)i * delta_beta);
// 		Weights.push_back(1.0);
//
// 		// assigne the number of steps counter
// 		state_visited.push_back(false);
// 		ave_V.push_back(new AverageObservable);
// 	}
//
// 	// initialise random number generator
// 	long unsigned seed;
// 	struct timeval tv;
// 	gettimeofday(&tv,0);
// 	seed = tv.tv_sec + tv.tv_usec + 7 * omp_get_thread_num();
//
// 	uniform = gsl_rng_alloc(gsl_rng_mt19937);
// 	gsl_rng_set(uniform, seed);
// }
//
// void SimulatedTempering::normalise_weights()
// {
// 	// dereference the constant
//     double weight_sum = 0.0;
//
// 	// sum over the temperatures
//     for(unsigned i=0; i<number_of_temperatures; i++)
// 		weight_sum += Weights[i];
//
//     // set the rescaled weights
//     for(unsigned i=0; i<number_of_temperatures; i++)
// 		Weights[i] /= weight_sum;
// }
//
// void SimulatedTempering::update_temperature(Molecule* molecule_pt,
// 											const unsigned& step)
// {
// 	// if we're on a switch step then swith
// 	// otherwise do nothing
// 	if(step % switch_frequency == 0)
// 	{
// 		// make proposed index
// 		double V = *molecule_pt->potential_pt();
// 		double U = 0.0;
// 		int shift = 0;
//
// 		// check if the current step has been visited
// 		if(state_visited[current_index] == false)
// 			state_visited[current_index] = true;
//
// 		// observe the potential energy for this state
// 		ave_V[current_index]->observe(V);
//
// 		// calculate the weights for this step
// 		for(unsigned i=1; i<number_of_temperatures; i++)
// 		{
// 			// If a state has not been visited, assume it has the same PE as preceding state
//             // to increase transition probability (see http://dx.doi.org/10.1021/acs.jpcb.5b03381 )
// 			if(state_visited[i] != false)
// 				Weights[i] = Weights[i-1] * exp(0.5 * (beta[i] - beta[i-1])
// 					   * (ave_V[i]->get_average() + ave_V[i-1]->get_average()));
// 			else
// 				Weights[i] = Weights[i-1] * exp(0.5 * (beta[i] - beta[i-1])
// 												* ave_V[i-1]->get_average());
// 		}
//
// 		// generator uniform random variable
// 		U = gsl_rng_uniform(uniform);
//
// 		// propose shift up
// 		if(U >= 0.5 && current_index + 1 < number_of_temperatures)
// 			shift = 1;
//
// 		// propose shift down
// 		if(U < 0.5 && current_index > 0)
// 			shift = -1;
//
// 		// if the shift is non-zero
// 		if(shift != 0)
// 		{
// 			// calculate acceptance probability
// 			double accept = Weights[current_index + shift]
// 							/ Weights[current_index]
// 				* exp(-(beta[current_index + shift] - beta[current_index]) * V);
//
// 			 // if we accept the state, update the temperatures
// 		  	// and the index - else do nothing
// 		  	if(U <= accept)
// 		  	{
// 		  		// update temperature index
// 		  		current_index += shift;
//
// 				// rescale the momentum
// 				rescale_momentum(molecule_pt, shift);
//
// 		  		// update temperature
// 		  		molecule_pt->set_beta(beta[current_index]);
// 		  	}
// 		}
// 	}
// }
//
// void SimulatedTempering::rescale_momentum(Molecule* molecule_pt, const int& shift)
// {
// 	// get the previous index
// 	int previous_index = current_index - shift;
// 	unsigned number_of_particles = molecule_pt->nparticle();
// 	unsigned dim = molecule_pt->dim();
// 	double* p = NULL;
//
// 	// calculate the scale
// 	double scale = sqrt(beta[previous_index] / beta[current_index]);
//
// 	Particle* particle_i = NULL;
//
// 	for(unsigned i=0; i<number_of_particles; i++)
// 		for(unsigned j=0; j<dim; j++)
// 		{
// 			particle_i = molecule_pt->particle_pt(i);
// 			p = particle_i->p_pt(j);
//
// 			(*p) *= scale;
// 		}
// }
