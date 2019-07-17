#include "Molecules.hpp"

/******************************************************************************
                           Generic Molecule Class
 *****************************************************************************/
Molecule::~Molecule()
{
    for(unsigned i=0; i<number_of_particles; i++)
        delete Particles.at(i);
    // Particles.clear();
}

// return dimension
unsigned Molecule::dim() {return DIM;}

// return number of particles in molecule
unsigned Molecule::nparticle() {return number_of_particles;}

// return temperature
double Molecule::kt() {return kT;}

// return inverse temperature
double Molecule::beta() {return Beta;}

// return pointer to the i-th particle
Particle& Molecule::particle(const unsigned& i) {return *Particles.at(i);}

// return pointer to the potential
double& Molecule::potential() {return V;}

// return the potential
double Molecule::potential() const {return V;}

// return pointer to the laplacian
double& Molecule::laplace() {return L;}

// return the laplacian
double Molecule::laplace() const {return L;}

// set temperature
void Molecule::set_temperature(const double& temperature)
{
    kT = temperature;
    Beta = 1.0 / temperature;
}

// set the beta
void Molecule::set_beta(const double& beta)
{
    kT = 1.0 / beta;
    Beta = beta;
}

/******************************************************************************
                           Singelton Molecule Class
 *****************************************************************************/
// constructor
Singelton::Singelton(const Vector& q_0, const Vector& p_0, const Matrix& m,
                     const double& kt, const unsigned& dim) :
                        Molecule(kt, 1, dim)
{
    // create one instance of a particle
    Particles.insert(make_pair(0, new Particle(q_0, p_0, m, dim)));
    // Particles[0] = new Particle(q_0, p_0, m, dim);
}

/******************************************************************************
                           Collection Molecule Class
 *****************************************************************************/
// class for a cluster of non-interactring particles
Collection::Collection(const Vector& q_0, const Vector& p_0, const Matrix& m,
                       const double& kt, const unsigned& nparts,
                       const unsigned& dim) :
                         Molecule(kt, nparts, dim)
{
    // create all the particles instance of a particle
    for(unsigned i=0; i<number_of_particles; i++)
        Particles.insert(make_pair(i, new Particle(q_0, p_0, m, dim)));
        // Particles[i] = new Particle(q_0, p_0, m, dim);
}

/******************************************************************************
                   Anisotropic Collection Molecule Class
 *****************************************************************************/
 AniCollection::AniCollection(const Vector& q_0, const Vector& p_0,
                              const Matrix& Q_0, const Matrix& pi_0,
                              const Matrix& m_0, const Matrix& I_0,
                              const unsigned& narms, const double& kt,
                              const unsigned& nparts, const unsigned& dim) :
                                Molecule(kt, nparts, dim)
{
    // create all the particles instance of a particle
    for(unsigned i=0; i<number_of_particles; i++)
        Particles.insert(make_pair(i, new Particle(q_0, p_0, Q_0, pi_0, m_0,
                                                   I_0, narms, dim)));
        // Particles[i] = new Particle(q_0, p_0, Q_0, pi_0, m_0, I_0,
                                      // narms, dim);
}

//
// // return the pointer to the laplacian of the potential
// double* Molecule::laplacian_potential_pt(){return &nabla_two_V;}
//
// // initialise the ISST beta weights
// void Molecule::ISST_set_beta_weight(const double& weight)
// {
//   // initialise the beta weight
//   ISST_Beta_Weight.push_back(weight);
//
//   // initialise the inverse beta weights
//   ISST_Beta_Inverted_Weight.push_back(1.0/weight);
//
//   // initialise the momentum
//   ISST_Beta_Weight_Momentum.push_back(0.0);
//
//   // initialise the force
//   ISST_Beta_Weight_Force.push_back(0.0);
//
//   // initalise the average sums
//   ISST_AverageSums.push_back(0.0);
// };
//
// // return the beta_weight vector
// double* Molecule::ISST_beta_w_pt(const double& i) {return &ISST_Beta_Weight[i];}
//
// // return pointer for Beta Weight momentum
// double* Molecule::ISST_beta_wp_pt(const double& i) {return &ISST_Beta_Weight_Momentum[i];}
//
// // return pointer for Beta Weight force
// double* Molecule::ISST_beta_w_force_pt(const double& i) {return &ISST_Beta_Weight_Force[i];}
//
// // return pointer for inverse beta weight
// double* Molecule::ISST_beta_inv_w_pt(const double& i) {return &ISST_Beta_Inverted_Weight[i];}
//
// // return the pointer to the i^th sum and j^th value
// double* Molecule::ISST_get_average_sum_pt(const unsigned& i)
// {
//   if (i > ISST_AverageSums.size())
//     {
//       printf("ERROR: ISST_AverageSums only holds %ld sums, trying to access %d",ISST_AverageSums.size(),i);
//       exit(-1);
//     }
//
//   // return the sum
//   return &ISST_AverageSums[i];
// }
//
// // set label j of particle i to LABEL
// void Molecule::set_particle_label(const unsigned& i,
// 				  const unsigned& j,
// 				  const int& LABEL)
// {
//   // set the j^th label of the i^th particle to
//   // be LABEL
//   Particle_pt[i]->set_label(j,LABEL);
// }
//
// // get label j of partile i
// int Molecule::get_particle_label(const unsigned& i,
// 				 const unsigned& j)
// {
//   // return the j^th label of the i^th particle
//   return Particle_pt[i]->get_label(j);
// }
//
// // set the posibility to do interpolation with
// // the molecule
// void Molecule::initialise_interpolation(const unsigned& nPnts)
// {
//   // set the number of interpolation points to two
//   for(unsigned i=0;i<Nparticle;i++)
//     Particle_pt[i]->init_interpolation(nPnts);
// };
//
// // update the interpolation of particle i
// void Molecule::update_particle_interpolation(const unsigned& i,
// 					     const double& q_n,
// 					     const double& p_n)
// {
//   //update interpolation on particle i
//   Particle_pt[i]->update_interpolation(q_n,p_n);
// };
//
// // return the weight
// double Molecule::get_W_T() {return W_T;}
//
// /******************************************************************************
//                            Singelton Molecule Class
//  *****************************************************************************/
// // constructor
// Singelton::Singelton(const std::vector<double>& q_0,
// 		     const std::vector<double>& p_0,
// 		     const std::vector<double>& M,
// 		     const double& kt,
// 		     const unsigned& dim)
// {
//   // set the number of particles to 1, dimension,
//   // inverse temperature,
//   Nparticle = 1;
//   Dim       = dim;
//
//   kT        = kt;
//   Beta      = 1.0/kt;
//
//   // create one instance of a particle
//   Particle_pt.resize(1,new Particle(dim));
//
//   // set global index
//   Particle_pt[0]->set_global_index(0);
//
//   // pass the initial conditions and the mass given to the particle
//   for(unsigned i=0;i<Dim;i++)
//     {
//       *Particle_pt[0]->q_pt(i) = q_0[i];
//       *Particle_pt[0]->p_pt(i) = p_0[i];
//       *Particle_pt[0]->m_pt(i) = M[i];
//     }
//
//   // set the holder for the weight in Paper with B-J-A-E
//   Molecule::W_T = 0.0;
// }
//
//
//
// /******************************************************************************
//                            Collection Molecule Class
//  *****************************************************************************/
// // class for a cluster of non-interactring particles
// Collection::Collection(const std::vector<double>& q_0,
// 		       const std::vector<double>& p_0,
// 		       const std::vector<double>& M,
// 		       const unsigned& nparticles,
// 		       const double& kt)
//
// {
//   // set basic parameters
//   Nparticle = nparticles;
//   Dim        = q_0.size();
//
//   kT         = kt;
//   Beta       = 1.0/kt;
//
//   // create all the particles instance of a particle
//   for(unsigned i=0;i<Nparticle;i++)
//     {
//       // make new particle
//       Particle_pt.push_back(new Particle(Dim));
//
//       // set 2 particle lables based on initiaal
//       // condition given in q_0[0] direction
//       Particle_pt[i]->initialise_label(2,(int)q_0[0]);
//
//       // pass the initial conditions and the mass given to the particle
//       for(unsigned j=0;j<Dim;j++)
// 	{
// 	  *Particle_pt[i]->q_pt(j) = q_0[j];
// 	  *Particle_pt[i]->p_pt(j) = p_0[j];
// 	  *Particle_pt[i]->m_pt(j) = M[j];
// 	}
//     }
// };
//
// // function which can be called to initialise the collection
// // in parallel for speed up.
// void Collection::initialise_parallel_performance(const std::vector<double>& q_0,
// 						 const std::vector<double>& p_0,
// 						 const std::vector<double>& M,
// 						 const unsigned& nparticle)
// {
//   // re-create all the particles on the individual
//   // processes
// #pragma omp for schedule(static)
//   for(unsigned i=0;i<Nparticle;i++)
//     {
//       // make new particle
//       Particle_pt[i] = new Particle(Dim);
//
//       // set 2 particle lables based on initiaal
//       // condition given in q_0[0] direction
//       Particle_pt[i]->initialise_label(2,(int)q_0[0]);
//
//       // pass the initial conditions and the mass given to the particle
//       for(unsigned j=0;j<Dim;j++)
// 	{
// 	  *Particle_pt[i]->q_pt(j) = q_0[j];
// 	  *Particle_pt[i]->p_pt(j) = p_0[j];
// 	  *Particle_pt[i]->m_pt(j) = M[j];
// 	}
//     }
//
// };
//
//
//
//
//
// /******************************************************************************
//                           Crystal Molecule Class
//  *****************************************************************************/
// // constructor
// Crystal::Crystal(const std::vector<double>& q_0,
// 		 const std::vector<double>& p_0,
// 		 const std::vector<double>& M,
// 		 const unsigned& nparticles,
// 		 const double& kt)
// {
//   // set the number of particles to 1, dimension,
//   // inverse temperature,
//   Nparticle = nparticles;
//   Dim        = q_0.size();
//
//   kT         = kt;
//   Beta       = 1.0/kt;
//
//   // create all the particles instance of a particle
//   for(unsigned i=0;i<Nparticle;i++)
//     {
//       // make new particle
//       Particle_pt.push_back(new Particle(Dim));
//
//       // set global index
//       Particle_pt[i]->set_global_index(i);
//
//       // pass the initial conditions and the mass given to the particle
//       for(unsigned j=0;j<Dim;j++)
// 	{
// 	  *Particle_pt[i]->q_pt(j) = q_0[j];
// 	  *Particle_pt[i]->p_pt(j) = p_0[j];
// 	  *Particle_pt[i]->m_pt(j) = M[j];
// 	}
//
//     }
//
//   // set the holder for the weight in Paper with B-J-A-E
//   Molecule::W_T = 0.0;
// }
//
// void Molecule::SetNeighbours()
// {
//   // if more neighbours we need more exclusions or figure out how to do it properly
//   for(unsigned i=0;i<Nparticle;i++)
//     {
//       if(i == Nparticle-2)
//        	{
//        	  (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[Nparticle-1]);
//       	  (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[0]);
//       	}
//        else if(i == Nparticle-1)
//       	 {
// 	   (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[0]);
//       	   (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[1]);
//       	 }
//        else
// 	 {
// 	   (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[i+1]);
// 	   (*Particle_pt[i]->neighbour_list()).push_back(Particle_pt[i+2]);
// 	 }
//     }
// }
//
// /******************************************************************************
//                    Anisotropic Crystal Molecule Class
//  *****************************************************************************/
// // constructor
// AniCrystal::AniCrystal(const std::vector<double>& q_0,
// 		       const std::vector<double>& p_0,
// 		       const std::vector<double>& M,
// 		       const double& pi,
// 		       const double& I,
// 		       const unsigned& nArms,
// 		       const unsigned& nparticles,
// 		       const double& kt)
// {
//   // set the number of particles to 1, dimension,
//   // inverse temperature,
//   Nparticle = nparticles;
//   Dim        = q_0.size();
//
//   kT         = kt;
//   Beta       = 1.0/kt;
//
//   // create all the particles instance of a particle
//   for(unsigned i=0;i<Nparticle;i++)
//     {
//       Particle_pt.push_back(new Particle(Dim));
//
//       // set global index
//       Particle_pt[i]->set_global_index(i);
//
//       // needs to be called before setting I and pi
//       Particle_pt[i]->set_eqidistant_arms(nArms);
//
//       *Particle_pt[i]->I_pt(0)  = I;
//       *Particle_pt[i]->pi_pt(0) = pi;
//
//       // pass the initial conditions and the mass given to the particle
//       for(unsigned j=0;j<Dim;j++)
// 	{
// 	  *Particle_pt[i]->q_pt(j) = q_0[j];
// 	  *Particle_pt[i]->p_pt(j) = p_0[j];
// 	  *Particle_pt[i]->m_pt(j) = M[j];
// 	}
//     }
//
//   // set the holder for the weight in Paper with B-J-A-E
//   Molecule::W_T = 0.0;
//  }
