#include "Particle.hpp"

using namespace::std;

/******************************************************************************
                           Generic Particle Class
 *****************************************************************************/
Particle::~Particle()
{
    // for(unsigned i=0; i<DIM; i++)
    // {
    //     delete[] &q[i];
    //     delete[] &p[i];
    //     delete[] &f[i];
    //     delete[] &m[i];
    // }

    // q.clear();
    // p.clear();
    // f.clear();
    // m.clear();
}

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
void Particle::set_eqidistant_arms()
{
    double narms = (double)arm.size()[1];
    double TwoPi = 2 * M_PI;

    for(unsigned i=0; i<arm.size()[1]; i++)
        if(DIM == 2)
        {
            arm(i, 0) = cos((double)i * TwoPi / narms);
            arm(i, 1) = sin((double)i * TwoPi / narms);
        }
}
