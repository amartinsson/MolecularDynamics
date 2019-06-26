#include "Particle.hpp"

using namespace::std;

/******************************************************************************
                           Generic Particle Class
 *****************************************************************************/
Particle::~Particle() {}

// return rigid body state
bool Particle::rigid_body() const
{
    return rigid_body_state;
}

// access to the dimension of the particle
unsigned Particle::dim() const
{
    return DIM;
}

// set equidistant arms
void Particle::set_eqidistant_arms(const unsigned& N)
{
    double narms = (double)N;
    double TwoPi = 2 * M_PI;

    for(unsigned i=0; i<N; i++)
    {
        arm(i, 0) = cos((double)i * TwoPi / narms);
        arm(i, 1) = sin((double)i * TwoPi / narms);
    }
}
