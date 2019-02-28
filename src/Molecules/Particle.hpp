#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <iostream>
#include <math.h>
#include <vector>

#include "Matrix.hpp"

using namespace::std;


/******************************************************************************
                           Generic Particle Class
 *****************************************************************************/
class Particle
{
public:
    // constructor
    Particle(const unsigned& dim);

    // destructor
    ~Particle();

    // access pointer to position index i
    double* q_pt(const unsigned i);

    // access pointer to momentum to position index i
    double* p_pt(const unsigned i);

    // access pointer to force
    double* f_pt(const unsigned i);

    // access pointer to mass
    double* m_pt(const unsigned i);

    // access to the pointer of the i,j element of the rotation matrix
    //double* Q_pt(const unsigned& i, const unsigned& j);
    double Q(const int& i, const int& j) const;
    double& Q(const int& i, const int& j);

    Matrix Q() const;
    Matrix& Q();

    // access to the pointer of the ith arm
    double arm(const unsigned& i, const unsigned& j);

    // access pointer angular momentum of interia index i
    double* pi_pt(const unsigned i);

    // access pointer to moment of interia index i
    double* I_pt(const unsigned i);

    // access pointer to torque index i
    double* tau_pt(const unsigned i);

    // return the state of the rigid body
    bool rigid_body();

    // return the dimesnion of the particle
    unsigned dim();

    // set equidistant arms
    void set_eqidistant_arms(const unsigned& narms);

private:
    // position
    std::vector<double> q;
    // momentum
    std::vector<double> p;
    // force
    std::vector<double> f;
    // mass
    std::vector<double> m;
    // rotation matrix
    Matrix Qmat;
    // vector of arms
    std::vector<double> arms;
    // anglular momentum
    std::vector<double> pi;
    // moment of intertia
    std::vector<double> I;
    // torque
    std::vector<double> tau;
    // holder of dimension
    unsigned DIM;
    // holder for the rigid body state
    bool rigid_body_state;
};

#endif
