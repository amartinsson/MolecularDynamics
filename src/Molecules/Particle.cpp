#include "Particle.hpp"

using namespace::std;

/******************************************************************************
                           Generic Particle Class
 *****************************************************************************/
Particle::Particle(const unsigned& dim) : Qmat(2, 2)
{
    // initialise all the vectors
    q.resize(dim, 0.0);
    p.resize(dim, 0.0);
    f.resize(dim, 0.0);
    m.resize(dim, 0.0);

    // set the dimension
    DIM = dim;

    // set the default rigid body state to zero
    rigid_body_state = false;
}

Particle::~Particle()
{
    for(unsigned i=0; i<DIM; i++)
    {
        delete[] &q[i];
        delete[] &p[i];
        delete[] &f[i];
        delete[] &m[i];
    }

    q.clear();
    p.clear();
    f.clear();
    m.clear();
}

// access pointer to position index i
double* Particle::q_pt(const unsigned i)
{
  return &q[i];
}

// access pointer to momentum index i
double* Particle::p_pt(const unsigned i)
{
  return &p[i];
}

// access pointer to force index i
double* Particle::f_pt(const unsigned i)
{
  return &f[i];
}

// access pointer to mass index i
double* Particle::m_pt(const unsigned i)
{
  return &m[i];
}

// access pointer to the i,j element of the rotation matrix
double Particle::Q(const int& i, const int& j) const
{
    return Qmat(i, j);
}

// access for updating the matrix
double& Particle::Q(const int& i, const int& j)
{
    return Qmat(i,j);
}

// access to the matrix
Matrix Particle::Q() const
{
    return Qmat;
}

// access to the entire rotation matrix
Matrix& Particle::Q()
{
    return Qmat;
}

// access to the jth component of the ith arm
double Particle::arm(const unsigned& j, const unsigned& i)
{
  return arms[j + i * DIM];
}

// access pointer angular momentum index i
double* Particle::pi_pt(const unsigned i)
{
  return &pi[i];
}

// access pointer to moment of interia index i
double* Particle::I_pt(const unsigned i)
{
  return &I[i];
}

// access pointer to torque index i
double* Particle::tau_pt(const unsigned i)
{
  return &tau[i];
}

// return rigid body state
bool Particle::rigid_body()
{
    return rigid_body_state;
}

// access to the dimension of the particle
unsigned Particle::dim()
{
    return DIM;
}

// set equidistant arms
void Particle::set_eqidistant_arms(const unsigned& narms)
{
  if(DIM != 2)
    {
        std::cout << "ERROR: CANNOT USE EQUIDISTANT ARMS WITH MORE "
                  << "THAN 2 dimensions" << std::endl;
        exit(-1);
    }

    // set the state to incorporate rigid body dynamics
    rigid_body_state = true;

    // resize the vectors and matrices
    // Q = Matrix(2,2);
    // Q.resize(DIM*DIM, 0.0);
    arms.resize(3*DIM, 0.0); // senare only 3 arms in particles

    if(DIM == 2)
    {
        I.resize(1,0.0);
        pi.resize(1,0.0);

        tau.resize(1,0.0);
    }
    else if(DIM == 3)
    {
        I.resize(DIM,0.0);
        pi.resize(DIM,0.0);

        tau.resize(DIM,0.0);
    }

    double TwoPi = 2 * M_PI;

    for(unsigned i=0;i<narms;i++)
    {
        arms[0 + i * DIM] = cos((double)i*TwoPi/(double)narms);
        arms[1 + i * DIM] = sin((double)i*TwoPi/(double)narms);
    }
}
