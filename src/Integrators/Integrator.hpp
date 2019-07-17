#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

#include "Molecules.hpp"

using namespace::std;

 /******************************************************************************
                            Generic Base Class
  *****************************************************************************/
class Integrator
{
public:
    Integrator() {};
    // destructor
    ~Integrator() {};
    // must be implemented
    virtual void integrate(Molecule* molecule_pt) = 0;
};

#endif
