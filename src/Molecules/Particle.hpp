#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <iostream>
#include <math.h>
#include <vector>

#include "Array.hpp"

using namespace::std;


/******************************************************************************
                           Generic Particle Class
 *****************************************************************************/
class Particle
{
public:
    // empty constructor
    Particle() : DIM(0), rigid_body_state(false) {};

    // explicit constructor
    explicit Particle(const unsigned& dim) : q(dim), p(dim), f(dim),
                                             m(dim, dim),
                                             rigid_body_state(false),
                                             DIM(dim) {}
    // explicit constructor
    explicit Particle(const Vector& q_0, const Vector& p_0, const Matrix& m_0,
                      const unsigned& dim) : q(q_0), p(p_0), f(dim), m(m_0),
                                             rigid_body_state(false),
                                             DIM(dim) {}
    // explicit constructor
    explicit Particle(const Vector& q_0, const Vector& p_0,
                      const Matrix& Q_0, const Matrix& pi_0,
                      const Matrix& m_0, const Matrix& I_0,
                      const unsigned& narms, const unsigned& dim) :
                      q(q_0), p(p_0), f(dim), Q(Q_0), pi(pi_0),
                      tau(1, 1),
                      m(m_0), I(I_0), rigid_body_state(true),
                      arm(narms, dim), DIM(dim)
    {
        // set the arms of the particle to correct number
        this->set_eqidistant_arms(narms);
    }
    // destructor
    ~Particle();

    // particle data
    Vector q;
    Vector p;
    Vector f;
    Matrix m;
    double laplace;

    // rotation data
    Matrix Q;
    Matrix pi;
    Matrix tau;
    Matrix I;

    // rotation arms
    Matrix arm;

    void operator= (const Particle& particle)
    {
        this->q = particle.q;
        this->p = particle.p;
        this->f = particle.f;
        this->m = particle.m;

        this->Q = particle.Q;
        this->pi = particle.pi;
        this->tau = particle.tau;
        this->I = particle.I;
        this->arm = particle.arm;

        this->rigid_body_state = particle.rigid_body();
        this->DIM = particle.dim();
    }

    // return the state of the rigid body
    bool rigid_body() const;

    // return the dimesnion of the particle
    unsigned dim() const;

private:
    // holder of dimension
    unsigned DIM;
    // holder for the rigid body state
    bool rigid_body_state;

    // set equidistant arms
    void set_eqidistant_arms(const unsigned& i);
};

#endif
