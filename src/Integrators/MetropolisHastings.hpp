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
    // potential
    double V;
    double Vnp1;
    // normal generator to generate random numbers
    NormalGenerator normal_gen;
    // uniform generator
    UniformGenerator uniform_gen;
    // system holder
    System* system;
    // boolean for check if first step
    bool first_step;

    // proposal for the algorithm
    Vector proposal_move(const Vector& q);
    // proposal ratio
    double proposal_ratio(Molecule* molecule_pt);
    // make an accept reject evaulation with acceptance a
    bool accept_reject(const double& a, const double& proposal_ratio);
    // creat the map
    void create_position_map(Molecule* molecule_pt);

private:
    // private
    AverageObservable* acceptance;
    // system for storing positions
    unordered_map<Particle*, Vector> qprime;
};

#endif
