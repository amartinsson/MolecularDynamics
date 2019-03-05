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
    Molecule(const double& kt, const unsigned& nparts, const unsigned& dim)
        : DIM(dim), number_of_particles(nparts), kT(kt), Beta(1.0 / kt),
          V(0.0) {Particles.resize(number_of_particles);}

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
    Particle& particle(const unsigned& i);
    // return pointer to the potential
    double& potential();
    // set temperature
    void set_temperature(const double& temperature);
    // set the beta
    void set_beta(const double& beta);
    // vector of pointers to particles
    vector<Particle*> Particles;

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
};

/******************************************************************************
                           Singelton Molecule Class
 *****************************************************************************/
class Singelton : public Molecule
{
public:
    Singelton(const Vector& q_0, const Vector& p_0, const Matrix& m,
              const double& kt, const unsigned& dim);
};


/******************************************************************************
                            Collection Molecule Class
 *****************************************************************************/
class Collection : public Molecule
{
public:
    Collection(const Vector& q_0, const Vector& p_0, const Matrix& m,
               const double& kt, const unsigned& nparts, const unsigned& dim);
};

/******************************************************************************
                   Anisotropic Collection Molecule Class
 *****************************************************************************/
class AniCollection : public Molecule
{
public:
    AniCollection(const Vector& q_0, const Vector& p_0,
                  const Matrix& Q_0, const Matrix& pi_0,
                  const Matrix& m_0, const Matrix& I_0,
                  const unsigned& narms, const double& kt,
                  const unsigned& nparts, const unsigned& dim);
};

#endif
