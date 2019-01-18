#ifndef HARMONICOSCILLATOR_HPP
#define HARMONICOSCILLATOR_HPP

#include "System.hpp"

using namespace::std;

/******************************************************************************
                           Harmonic Oscillator
 *****************************************************************************/
class HarmonicOscillator : public System
{
public:
    HarmonicOscillator(Molecule* molecule_pt);
    // destructor
    ~HarmonicOscillator(){};
    // compute the force
    void compute_force(Molecule* molecule_pt);
    // compute the pair force
	vector<double> compute_pair_force(Molecule* molecule_pt,
									  Particle* particle_i,
				                      Particle* particle_j,
									  const double& r, const double& r_x,
                                      const double& r_y);
};

#endif
