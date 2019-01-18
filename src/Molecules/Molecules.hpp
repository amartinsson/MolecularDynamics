#ifndef MOLECULES_HPP
#define MOLECULES_HPP

#include <vector>

#include "Particle.hpp"

using namespace::std;

/******************************************************************************
                           Generic Molecule Class
 *****************************************************************************/
class Molecule
{
public:
    // constructor
    Molecule();

    // destructor
    ~Molecule();

    // return dimension
    unsigned dim();
    // return number of particles in molecule
    unsigned nparticle();
    // return temperature
    double kt();
    // return inverse temperature
    double beta();
    // return pointer to the i-th particle
    Particle* particle_pt(const unsigned& i);
    // return pointer to the potential
    double* potential_pt();
    // set temperature
    void set_temperature(const double& temperature);
    // set the beta
    void set_beta(const double& beta);

protected:
    // Dimension
    unsigned DIM;
    // number of particles
    unsigned number_of_particles;
    //temperatures
    double kT;
    double Beta;
    //potential
    double V;
    // vector of pointers to particles
    std::vector<Particle*> Particle_pt;
};

/******************************************************************************
                           Singelton Molecule Class
 *****************************************************************************/
class Singelton : public Molecule
{
public:
    Singelton(const std::vector<double>& q_0, const std::vector<double>& p_0,
	          const std::vector<double>& m, const double& kt,
	          const unsigned& dim);
};


/******************************************************************************
                            Collection Molecule Class
 *****************************************************************************/
class Collection : public Molecule
{
public:
    Collection(const std::vector<double>& q_0, const std::vector<double>& p_0,
	     const std::vector<double>& M, const unsigned& nparticles,
	     const double& kt);
};

/******************************************************************************
                          Crystal Molecule Class
 *****************************************************************************/
class Crystal : public Molecule
{
public:
    Crystal(const std::vector<double>& q_0, const std::vector<double>& p_0,
	        const std::vector<double>& M, const unsigned& nparticles,
	        const double& kt);
};

/******************************************************************************
                   Anisotropic Crystal Molecule Class
 *****************************************************************************/
class AniCrystal : public Molecule
{
public:
    AniCrystal(const std::vector<double>& q_0, const std::vector<double>& p_0,
	           const std::vector<double>& M, const double& pi, const double& I,
	           const unsigned& nArms, const unsigned& nparticles,
               const double& kt);
};

#endif
