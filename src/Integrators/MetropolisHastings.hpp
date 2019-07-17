#ifndef METROPOLISHASTINGS_HPP
#define METROPOLISHASTINGS_HPP

#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>

#include "Array.hpp"
#include "Generator.hpp"
#include "Integrator.hpp"
#include "System.hpp"
#include "AverageObservable.hpp"

using namespace::std;

 /******************************************************************************
                            Metropolis Hastings Base Class
  *****************************************************************************/
class MetropolisHastings : public Integrator
{
public:
    MetropolisHastings(const double& beta, const double& sigma,
        System* system_pt, const int& seed);
    // destructor
    ~MetropolisHastings();
    // must be implemented
    void integrate(Molecule* molecule_pt);

protected:
    // holder for inverse temperature
    double beta;
    // holder for the variance
    double sigma;
    // normal generator to generate random numbers
    NormalGenerator normal_gen;
    // uniform generator
    UniformGenerator uniform_gen;
    // system holder
    System* system;

    // proposal for the algorithm
    Vector proposal_move(const Vector& q);
    // proposal ratio
    double proposal_ratio(const Vector& q, const Vector& qprime);
    // calculate the acceptance ratio
    double acceptance_ratio(const Vector& q, const Vector& qprime);

private:
    AverageObservable* acceptance;
};

#endif
