#include "OngulatedTempering.hpp"

using namespace::std;

/******************************************************************************
                            Ongulated Tempering
 *****************************************************************************/
// constructor
OngulatedTempering::OngulatedTempering(const double& tmin, const double& tmax,
	const double& n_temperatures, const unsigned& mod_switch, const int& seed)
		: SimulatedTempering(tmin, tmax, n_temperatures, mod_switch, seed) {};

// update the temperature weight
void OngulatedTempering::update_temperature(Molecule* molecule_pt,
	const unsigned& step) {

		// update temperature index
		int new_index = floor((cos(2.0 * M_PI / (double)switch_frequency * (double)step) + 1.0) * ((double)SimulatedTempering::number_of_temperatures-1.0) * 0.5);

		// index magic
		int shift = new_index - current_index;
		current_index = new_index;

		if(shift == 1 or shift == -1) {
			printf("ni=%d, ci=%d, s=%d, on step %d\n", new_index, current_index, shift, step);
			// rescale the momentum
			SimulatedTempering::rescale_momentum(molecule_pt, shift);

			// update temperature
			molecule_pt->set_beta(beta[current_index]);
		}
}
