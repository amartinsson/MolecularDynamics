#include "System.hpp"

using namespace::std;



// /******************************************************************************
//                               Generic System
//  *****************************************************************************/
// // System Constructor
// Systems::Systems(){};
//
// // function for the Lennard Jones potential
// double Systems::LJ_potential(const double& r)
// {
//   return 4 * (*Epsilon_pt) * (pow((*Sigma_pt)/r,12)
// 			      -pow((*Sigma_pt)/r,6));
// }
//
// std::vector<double> Systems::LJ_force(Particle* PartI,
// 		       Particle* PartJ,
// 		       const double& r,
// 		       const double& r_x,
// 		       const double& r_y)
// {
//   // tracker of the pair force
//   std::vector<double> forces(2, 0.0);
//
//   // ------------ LENNARD JONES FORCE
//   //
//   //    F_x = - 24 * eps * (2 * sigma^12 * r^-13 - sigma^6 r^-7) * x/r;
//   //    F_y = - 24 * eps * (2 * sigma^12 * r^-13 - sigma^6 r^-7) * y/r;
//   //
//
//   double dphidr = -24.*(*Epsilon_pt)*pow((*Sigma_pt),6.0)*
//     (2*pow((*Sigma_pt),6.0)*pow(r,-14.0) - pow(r,-8.0)) * r;
//
// // #pragma omp atomic
// //   *PartI->f_pt(0) += 2.0 * (dphidr * r_x / r);
// //
// // #pragma omp atomic
// //   *PartI->f_pt(1) += 2.0 * (dphidr * r_y / r);
// //
// // #pragma omp atomic
// //   *PartJ->f_pt(0) -= 2.0 * (dphidr * r_x / r);
// //
// // #pragma omp atomic
// //   *PartJ->f_pt(1) -= 2.0 * (dphidr * r_y / r);
//
// #pragma omp atomic
//   *PartI->f_pt(0) += dphidr * r_x / r;
//
// #pragma omp atomic
//   *PartI->f_pt(1) += dphidr * r_y / r;
//
// #pragma omp atomic
//   *PartJ->f_pt(0) -= dphidr * r_x / r;
//
// #pragma omp atomic
//   *PartJ->f_pt(1) -= dphidr * r_y / r;
//
//   // add and return the pair force
//   forces[0] = dphidr * r_x / r;
//   forces[1] = dphidr * r_y / r;
//
//   return forces;
// }
//
// // initialise the system
// void Systems::system_initialise(std::vector<Molecule*>& Mol_pt,
// 				const unsigned& nCopies)
// {
//   std::cout << "ERROR: system_initialise function must be declared in System"
// 	    << std::endl;
//   exit(-1);
// }
//
// // notify when not overridden
// void Systems::compute_force(Molecule* Mol_pt)
// {
//   std::cout << "ERROR: void function compute_force(Molecule* Mol_pt), "
// 	    << "must be implemented in System!" << std::endl;
//   exit(-1);
// }
//
// // function for computing the potential
// void Systems::compute_potential(Molecule* Mol_pt)
// {
//   std::cout << "ERROR: void function compute_potential(Molecule* Mol_pt), "
// 	    << "must be implemented in System!" << std::endl;
//   exit(-1);
// }
//
// // function for determing the region
// void Systems::determine_region(std::vector<Molecule*>& Mol_pt)
// {
//   std::cout << "ERROR: void function determine_region"
// 	    << "(std::vector<Molecule*>& Mol_pt), "
// 	    << "is currently only implemeted for DoubleWell System!"
// 	    << std::endl;
//   exit(-1);
// }
//
// /******************************************************************************
//                            An-Harmonic Oscillator
//  *****************************************************************************/
// // initialise static variables
// double AnHarmonicOscillator::Eta = 0.0;
// double AnHarmonicOscillator::Epsilon = 0.0;
//
// // initialise static variables
// double AnHarmonicOscillator::gauss_mu    = 0.0;
// double AnHarmonicOscillator::gauss_sigma = 1.0;
//
// // system initialise
// void AnHarmonicOscillator::system_initialise(Molecule* MolPt,
// 					     const double& epsilon,
// 					     const double& eta)
// {
//   // set the epsilon variable
//   AnHarmonicOscillator::Epsilon = epsilon;
//
//   // set the eta variable
//   AnHarmonicOscillator::Eta = eta;
//
//   // initalise force and momentum
//   compute_force(MolPt);
//   compute_potential(MolPt);
// };
//
// // compute force with molecule
// void AnHarmonicOscillator::compute_force(Molecule* MolPt)
// {
//   // dereference
//   double* f;
//   double q = 0.0;
//
// #pragma omp parallel for schedule(static)
//   for(unsigned j=0;j<MolPt->nparticle();j++)
//     for(unsigned i=0;i<MolPt->dim();i++)
//       {
// 	f = MolPt->particle_pt(j)->f_pt(i);
// 	q = *MolPt->particle_pt(j)->q_pt(i);
//
// 	*f =-(q + pow(q,2.0) * (Epsilon + Eta * q));
//       }
// };
//
// // compute force with particle
// void AnHarmonicOscillator::compute_force(Particle* PartPt,
// 					 const unsigned& DIM)
// {
//   // dereference
//   double* f;
//   double q = 0.0;
//
// #pragma omp simd
//   for(unsigned i=0;i<DIM;i++)
//     {
//       f = PartPt->f_pt(i);
//       q = *PartPt->q_pt(i);
//
//       *f =-(q + pow(q,2.0) * (Epsilon + Eta * q));
//     }
// };
//
// // function for computing the potential
// void AnHarmonicOscillator::compute_potential(Molecule* MolPt)
// {
//   // reset the potential
//   *MolPt->potential_pt()=0.0;
//
//   // dereference
//   double q = 0;
//   double* V = MolPt->potential_pt();
//
// #pragma omp parallel for schedule(static)
//   for(unsigned j=0;j<MolPt->nparticle();j++)
//     for(unsigned i=0;i<MolPt->dim();i++)
//       {
// 	q = *MolPt->particle_pt(j)->q_pt(i);
//
// #pragma omp atomic
// 	*V += 1.0 / 2.0 * pow(q,2.0)
// 	  + 1.0 / 3.0 * Epsilon * pow(q,3.0)
// 	  + 1.0 / 4.0 * Eta * pow(q,4.0);
//       }
// };
//
// // function returning the partition function integrand at position x
// double AnHarmonicOscillator::z_beta(double x, void * params)
// {
//   // dereference
//   double beta = *(double *) params;
//
//   double V = 1.0 / 2.0 * pow(x,2.0)
//     + 1.0 / 3.0 * Epsilon * pow(x,3.0)
//     + 1.0 / 4.0 * Eta * pow(x,4.0);
//
//   // make function
//   double f = exp(-beta * V);
//   return f;
// }
//
// // function returning the partition function integrand at position x
// double AnHarmonicOscillator::q2_beta(double x, void * params)
// {
//   // dereference
//   double beta = *(double *) params;
//
//   double V = 1.0 / 2.0 * pow(x,2.0)
//     + 1.0 / 3.0 * Epsilon * pow(x,3.0)
//     + 1.0 / 4.0 * Eta * pow(x,4.0);
//
//   // make function
//   double f = x * x * exp(-beta * V);
//   return f;
// }
//
// double AnHarmonicOscillator::gauss_beta(double x, void * params)
// {
//   // dereference
//   double beta = *(double *) params;
//
//   double V = 1.0 / 2.0 * pow(x,2.0)
//     + 1.0 / 3.0 * Epsilon * pow(x,3.0)
//     + 1.0 / 4.0 * Eta * pow(x,4.0);
//
//   // make function
//   double f = exp(-0.5 * pow(x - gauss_mu, 2.0) * pow(gauss_sigma,-2.0)) * exp(-beta * V);
//   return f;
// }
//
//
// // return exact pdf
// double AnHarmonicOscillator::get_exact_pdf(const double& x,
// 					   const double& Beta)
// {
//   // dereference
//   double V     = 0.0;
//   double Z     = 0.0;
//   double beta  = Beta;
//
//   // calculate the potential at point q
//   V += 1.0 / 2.0 * pow(x,2.0)
//     + 1.0 / 3.0 * Epsilon * pow(x,3.0)
//     + 1.0 / 4.0 * Eta * pow(x,4.0);
//
//   // calculate the normalising constant
//   gsl_integration_workspace * w
//     = gsl_integration_workspace_alloc (10000);
//
//   double result, error;
//
//   gsl_function F;
//   F.function = &AnHarmonicOscillator::z_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//
//   Z = result;
//
//   // free the workspace
//   gsl_integration_workspace_free(w);
//
//   // return the pdf at point x
//   return exp(-beta * V)/ Z;
// };
//
//
// // return the exact solution to < q^2 >
// double AnHarmonicOscillator::get_exact_q2(const double& Beta)
// {
//   // dereference
//   double Z     = 0.0;
//   double beta  = Beta;
//
//   // -------------------------- NORMALISION CONSTANT
//   gsl_integration_workspace * w
//     = gsl_integration_workspace_alloc (10000);
//
//   double result, error;
//   gsl_function F;
//
//   F.function = &AnHarmonicOscillator::z_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//   Z = result;
//
//   // -------------------------- OBSERVBLE
//   double qTwo  = 0.0;
//
//   F.function = &AnHarmonicOscillator::q2_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//   qTwo = result;
//
//   // free the workspace
//   gsl_integration_workspace_free(w);
//
//   // return the pdf at point x
//   return qTwo / Z;
// };
//
// // return the exact solution to < exp(-0.5 * (q-mu)^2/sigma^2)) >
// double AnHarmonicOscillator::get_exact_gauss(const double& Beta,
// 					     const double& mu,
// 					     const double& sigma)
// {
//    // dereference
//   double Z     = 0.0;
//   double beta  = Beta;
//
//   AnHarmonicOscillator::gauss_mu    = mu;
//   AnHarmonicOscillator::gauss_sigma = sigma;
//
//   // -------------------------- NORMALISION CONSTANT
//   gsl_integration_workspace * w
//     = gsl_integration_workspace_alloc (10000);
//
//   double result, error;
//   gsl_function F;
//
//   F.function = &AnHarmonicOscillator::z_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//   Z = result;
//
//   // -------------------------- OBSERVBLE
//   double gauss  = 0.0;
//
//   F.function = &AnHarmonicOscillator::gauss_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//   gauss = result;
//
//   // free the workspace
//   gsl_integration_workspace_free(w);
//
//   // return the pdf at point x
//   return gauss / Z;
// };
//
// /******************************************************************************
//                                Double Well
//  *****************************************************************************/
// // constructor
// DoubleWell::DoubleWell(){};
//
// // initialise static doubles
// double DoubleWell::a = 0.0;
// double DoubleWell::b = 0.0;
//
// void DoubleWell::system_initialise(Molecule* Mol_pt,
// 				   double& A,double& B)
// {
//   // Check parameters and set parameters
//   if(B<A)
//     {
//       std::cout << "ERROR: Double Well parameters are b<a, "
// 		<< "with a=" << a << " and b=" << b << std::endl;
//       exit(-1);
//     }
//
//   // assign the parameters
//   DoubleWell::a = A;
//   DoubleWell::b = B;
//
//   // calculate initial force and set the initial markers
//   compute_force(Mol_pt);
// }
//
// // funciton for computing the force
// void DoubleWell::compute_force(Molecule* MolPt)
// {
//   // dereference
//   //  double* f;
//   double  q = 0.0;
//
//   unsigned Dim = MolPt->dim();
//   unsigned nParticles = MolPt->nparticle();
//
// #pragma omp for schedule(static)
//   for(unsigned i=0;i<nParticles;i++)
//     {
// #pragma omp simd
//       for(unsigned j=0;j<Dim;j++)
// 	{
// 	  double* f = MolPt->particle_pt(i)->f_pt(j);
// 	  q = *MolPt->particle_pt(i)->q_pt(j);
//
// 	  *f= -(4.0 * b - 2.0 * a) * ( q*q - 1.0) * q - 0.5 * a;
//       }
//     }
// }
//
// // function for computing the potential
// void DoubleWell::compute_potential(Molecule* MolPt)
// {
//   // dereference
//   double q  = 0.0;
//   double* V = MolPt->potential_pt();
//
//   unsigned Dim = MolPt->dim();
//   unsigned nParticles = MolPt->nparticle();
//
//   // reset the potentail
//   *V = 0.0;
//
//   // loop over all the dimensions
// #pragma omp simd collapse(2)
//   for(unsigned i=0;i<nParticles;i++)
//     for(unsigned j=0;j<Dim;j++)
//       {
// 	q = *MolPt->particle_pt(i)->q_pt(j);
//
// 	*V += (b - 0.5*a)*(q*q -1.0)*(q*q -1.0) + 0.5*a * (q + 1.0);
//       }
// }
//
// // function returning the partition function integrand at position x
// double DoubleWell::z_beta(double x, void * params)
// {
//   // dereference
//   double beta = *(double *) params;
//
//   double V = (b - 0.5*a)*(x*x -1.0)*(x*x -1.0)
//     + 0.5*a * (x + 1.0);
//
//
//   // make function
//   double f = exp(-beta * V);
//   return f;
// }
//
// // return exact pdf
// double DoubleWell::get_exact_pdf(const double& x,
// 				 const double& Beta)
// {
//   // dereference
//   double V     = 0.0;
//   double Z     = 0.0;
//   double beta  = Beta;
//
//   // calculate the potential at point q
//   V += (b - 0.5*a)*(x*x -1.0)*(x*x -1.0)
//     + 0.5*a * (x + 1.0);
//
//   // calculate the normalising constant
//   gsl_integration_workspace * w
//     = gsl_integration_workspace_alloc (10000);
//
//   double result, error;
//
//   gsl_function F;
//   F.function = &z_beta;
//   F.params = &beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//
//   Z = result;
//
//   // free the workspace
//   gsl_integration_workspace_free(w);
//
//   // return the pdf at point x
//   return exp(-beta * V)/ Z;
// };
//
// double DoubleWell::get_normalisation_constant(const double& beta)
// {
//   // dereference
//   double Z     = 0.0;
//   double Beta = beta;
//
//   // calculate the normalising constant
//   gsl_integration_workspace * w
//     = gsl_integration_workspace_alloc (1000);
//
//   double result, error;
//
//   gsl_function F;
//   F.function = &z_beta;
//   F.params = &Beta;
//
//   gsl_integration_qagi(&F,0,1e-7,1000,w, &result, &error);
//   Z = result;
//
//   // free the workspace
//   gsl_integration_workspace_free(w);
//
//   // return the pdf at point x
//   return Z;
// };
//
// /******************************************************************************
//                                Trimer
//  *****************************************************************************/
// // constructor
// Trimer::Trimer(){};
//
// // function for initilasing the system
// void Trimer::system_initialise(Molecule* MolPt,
// 			       const double& sigma,
// 			       const double& epsilon)
// {
//   // Initialise the sigma parameter for LJ-6,12
//   //Sigma_pt=&sigma;
//
//   // Initialise the epsilon parameter for LJ-6,12
//   //Epsilon_pt=&epsilon;
//
//   // initalise force and momentum
//   compute_force(MolPt);
//   compute_potential(MolPt);
// };
//
// // function for computing the force
// void Trimer::compute_force(Molecule* Mol_pt)
// {
//   // find the distance from the origin
//   // double r=sqrt((*Mol_pt->particle_pt(0)->q_pt(0))
//   //  		*(*Mol_pt->particle_pt(0)->q_pt(0))
//   //  		+(*Mol_pt->particle_pt(0)->q_pt(1))
//   //  		*(*Mol_pt->particle_pt(0)->q_pt(1)));
//
//   std::cout << "SENARE, CHANGE THE compute_force in Trimer used LJ_force previously" << std::endl;
//   exit(-1);
//   // compute the force in each direction
//   // *Mol_pt->particle_pt(0)->f_pt(0)=-2*(DL_DR)*(*Mol_pt->particle_pt(0)->q_pt(0))
//   //   -2*dLJ_dr(2*(*Mol_pt->particle_pt(0)->q_pt(0)));
//   // *Mol_pt->particle_pt(0)->f_pt(1)=-2*(DL_DR)*
//   //   (*Mol_pt->particle_pt(0)->q_pt(1));
// }
//
// /******************************************************************************
//                                   Vanden
//  *****************************************************************************/
// // Constructor
// Vanden::Vanden(){R=0.0;};
//
// // Function for initilasing the system
// void Vanden::system_initialise(Molecule* Mol_pt,
// 			       double& lambda)
// {
//   // Set pointer to the lambda strength vector
//   Lambda = lambda;
//
//   // Initialise the force and potential
//   compute_force(Mol_pt);
//   compute_potential(Mol_pt);
// };
//
// // Function returning the force scaling in Infinite Replica Exchange R_k
// double* Vanden::r_pt()
// {
//   return &R;
// }
//
// // Function for computing the force
// void Vanden::compute_force(Molecule* Mol_pt)
// {
//   // Helper parameters
//   unsigned DIM = Mol_pt->dim();
//   double* f;
//
//   double q = 0.0;
//
//   // Loop over all the dimensions, first dimension is special
//   for(unsigned i=0;i<DIM;i++)
//     {
//       // dereference particle
//       f = Mol_pt->particle_pt(0)->f_pt(i);
//       q = *Mol_pt->particle_pt(0)->q_pt(i);
//
//       if(i==0)
// 	{
// 	  *f = -4.0 * q * (1.0 - q*q) + 0.25;
// 	}
//       else
// 	{
// 	  *f = -Lambda * q;
// 	}
//     }
// };
//
// // function for computing the potential
// void Vanden::compute_potential(Molecule* Mol_pt)
// {
//   // Helper parameters
//   unsigned DIM = Mol_pt->dim();
//   double q = 0.0;
//
//   double V = 0.0;
//
//   // Loop over all the dimensions, first dimension is special
//   for(unsigned i=0;i<DIM;i++)
//     {
//       // dereference particle
//       q = *Mol_pt->particle_pt(0)->q_pt(i);
//
//       if(i==0)
// 	{
// 	  V += (1.0 - q*q) * (1.0 - q*q) - 0.25 * q;
// 	}
//       else
// 	{
// 	  V += 0.5 * Lambda * q*q;
// 	}
//     }
//
//   // update the potential
//   *Mol_pt->potential_pt() = V;
//
//   if(V != V)
//     {
//       printf("potential V = %f\n",V);exit(-1);
//     }
// };
//
// /******************************************************************************
//                                 Mercedes-Benz
//  *****************************************************************************/
// MercBenz::MercBenz(){};
//
// // special helper function for this potential
// double MercBenz::G(const double& x){return exp(-0.5*(x*x)/(Sigma_HB*Sigma_HB));};
//
// // Function for initialising the system
// void MercBenz::system_initialise(Molecule* MolPt,
// 				 double& sigma,
// 				 double& epsilon,
// 				 double& sigma_hb,
// 				 double& r_hb,
// 				 double& epsilon_hb,
// 				 const double& L_d)
// {
//   // check that the number of dimension is 2
//   if(2 != MolPt->dim())
//     {
//       printf("\nERROR: the dimension given to MercBenz is %i but should be 2\n\n",MolPt->dim());
//       exit(-1);
//     }
//
//   // set the sigma in the Lennard-Jones potential
//   Sigma_pt = &sigma;
//
//   // set the epsilon in the Lennard-Jones potential
//   Epsilon_pt = &epsilon;
//
//   // set sigma in the Hydrogen Bond part
//   Sigma_HB = sigma_hb;
//
//   // set radius Hydrogen Bond part
//   r_HB = r_hb;
//
//   // set epsilon Hydrogen Bond part
//   Epsilon_HB = epsilon_hb;
//
//   // make sure number of particles is perfect square
//   NumPartRow = sqrt(MolPt->nparticle());
//   NumPartCol = sqrt(MolPt->nparticle());
//
//   unsigned m = 0;
//   double SepX = 2.0*L_d/ (3.0*(double)NumPartRow);
//   //double SepX = L_d/((double)NumPartRow+ 1.0);
//
//   double SepY = L_d/((double)NumPartCol+ 1.0);
//   double phi_0 = 0.0;
//
//   double X = 0.0;
//   double Y = 0.0;
//
//   for(unsigned i=0;i<NumPartCol;i++)
//   {
//       for(unsigned j=0;j<NumPartRow;j++)
//       {
//           m = j + i * NumPartRow;
//
//        	  if(i % 2 == 0)
//           {
//               if(j==0)
//               {
//                   phi_0 = 0.0;
//                   X = 1.0 * SepX;
//                   Y = L_d - 0.5 * SepY - i * SepY;
//        		   }
//        	      else if(j % 2 ==0)
//               {
//                   phi_0 = 0.0;
//                   X += 2.0 * SepX;
//                   //X += 1.0 * SepX
//               }
//        	      else
//               {
//                   phi_0 = M_PI/3.0;
//                   X += 1.0 * SepX;
//        		  }
//        	  }
//       	   else
//            {
//                if( j == 0 )
//                {
//                    phi_0 = M_PI/3.0;
// 		           X = 0.5 * SepX;
// 		           Y = L_d - 0.5 * SepY - i * SepY;
//                }
//        	      else if(j % 2 ==0)
//        		  {
//                   phi_0 = M_PI/3.0;
// 		          X += 1.0 * SepX;
//        		  }
//        	      else
//        		  {
//                   phi_0 = 0.0;
// 		          X += 2.0 * SepX;
// 		          //X += 1.0 * SepX;
//        		  }
//   	       }
//
// 	  *MolPt->particle_pt(m)->q_pt(0) = X;
//   	  *MolPt->particle_pt(m)->q_pt(1) = Y;
//
//   	  *MolPt->particle_pt(m)->Q_pt(0,0) =  cos(phi_0);
//   	  *MolPt->particle_pt(m)->Q_pt(0,1) =  sin(phi_0);
//
//   	  *MolPt->particle_pt(m)->Q_pt(1,0) = -sin(phi_0);
//   	  *MolPt->particle_pt(m)->Q_pt(1,1) =  cos(phi_0);
//   	}
//     }
// };
//
// // function for computing the force
// std::vector<double> MercBenz::compute_pair_force(Molecule* MolPt,
//                                                  Particle* PartI,
// 				                                         Particle* PartJ,
// 				                                         const double& r,
// 				                                         const double& r_x,
// 				                                         const double& r_y)
// {
//   // holders for computing the forces
//   // between the particles
//   std::vector<double> lj_forces(2, 0.0);
//   std::vector<double> forces(2, 0.0);
//
//   // Lennard Jones Force and Potential
//   lj_forces = LJ_force(PartI, PartJ, r, r_x, r_y);
//
//   forces[0] += lj_forces[0];
//   forces[1] += lj_forces[1];
//
// #pragma omp atomic
//   *MolPt->potential_pt() += LJ_potential(r);
//
//   // ------------------------------ HYDROGEN BOND ----------------------------------- //
//
//   // Helpers and derivatives with respect to x and y
//   std::vector<double> h_i(3,0.0);
//   std::vector<double> h_j(3,0.0);
//
//   std::vector<double> dh_i_dx_i(3,0.0);
//   std::vector<double> dh_j_dx_i(3,0.0);
//
//   std::vector<double> dG_dh_i(3,0.0);
//   std::vector<double> dG_dh_j(3,0.0);
//
//   std::vector<double> dh_i_dy_i(3,0.0);
//   std::vector<double> dh_j_dy_i(3,0.0);
//
//   std::vector<double> dh_i_dQ_00(3,0.0);
//   std::vector<double> dh_i_dQ_01(3,0.0);
//
//   std::vector<double> dh_i_dQ_10(3,0.0);
//   std::vector<double> dh_i_dQ_11(3,0.0);
//
//   std::vector<double> dh_j_dQ_00(3,0.0);
//   std::vector<double> dh_j_dQ_01(3,0.0);
//
//   std::vector<double> dh_j_dQ_10(3,0.0);
//   std::vector<double> dh_j_dQ_11(3,0.0);
//
//   double dphiHB_dx = 0.0;
//   double dphiHB_dy = 0.0;
//
//   double dphiHB_dQ_i_00 = 0.0;
//   double dphiHB_dQ_i_01 = 0.0;
//
//   double dphiHB_dQ_i_10 = 0.0;
//   double dphiHB_dQ_i_11 = 0.0;
//
//   double dphiHB_dQ_j_00 = 0.0;
//   double dphiHB_dQ_j_01 = 0.0;
//
//   double dphiHB_dQ_j_10 = 0.0;
//   double dphiHB_dQ_j_11 = 0.0;
//
//   // derivatives with respect to r
//   double dG_dr = -(r-r_HB) / (Sigma_HB*Sigma_HB) * G(r-r_HB);
//   double dr_dx = -r_x/r;
//   double dr_dy = -r_y/r;
//
//   // unit vector in r direction
//   std::vector<double> u(2,0.0);
//   u[0] = r_x/r;
//   u[1] = r_y/r;
//
//   // local storage for particle
//   double Q_i[2][2];
//   double Q_j[2][2];
//
//   double arm_i[2][3];
//   double arm_j[2][3];
//
//   // read in orientation from memory
//   // for(unsigned k=0; k<2; k++)
//   //   for(unsigned l=0;l<2;l++)
//   //   {
//   //       // particle i
//   //       Q_i[k][l] = *PartI->Q_pt(k, l);
//   //
//   //       // particle j
//   //       Q_j[k][l] = *PartJ->Q_pt(k, l);
//   //   }
//   for(unsigned k=0; k<3; k++)
//   {
//         for(unsigned l=0;l<2;l++)
//         {
//             if(k<2)
//             {
//                 // particle i
//                 Q_i[k][l] = *PartI->Q_pt(k, l);
//
//                 // particle j
//                 Q_j[k][l] = *PartJ->Q_pt(k, l);
//             }
//
//             // read in particle i
//             arm_i[l][k] = PartI->arm(l, k);
//
//             // read in particle j
//             arm_j[l][k] = PartJ->arm(l, k);
//         }
//     }
//
//   // // read in arms from memory
//   // for(unsigned k=0; k<3; k++)
//   // {
//   //     // read in particle i
//   //     arm_i[0][k] = PartI->arm(0, k);
//   //     arm_i[1][k] = PartI->arm(1, k);
//   //
//   //     // read in particle j
//   //     arm_j[0][k] = PartJ->arm(0, k);
//   //     arm_j[1][k] = PartJ->arm(1, k);
//   // }
//
//   for(unsigned k=0; k<3; k++)
//   {
//       h_i[k] = (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) * u[0]
//              + (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) * u[1];
//
//       h_j[k] = (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) * u[0]
//            + (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) * u[1];
//
//       dh_i_dx_i[k] = r_x * h_i[k] / (r * r)
//                    - (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) / r;
//
//       dh_j_dx_i[k] = r_x * h_j[k] / (r * r)
//                    - (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) / r;
//
//       dh_i_dy_i[k] = r_y * h_i[k] / (r * r)
//                    - (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) / r;
//
//       dh_j_dy_i[k] = r_y * h_j[k] / (r * r)
//                    - (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) / r;
//
//       dh_i_dQ_00[k] = u[0] *  arm_i[0][k];
//       dh_i_dQ_01[k] = u[1] *  arm_i[0][k];
//
//       dh_i_dQ_10[k] = u[0] *  arm_i[1][k];
//       dh_i_dQ_11[k] = u[1] *  arm_i[1][k];
//
//        // matrix elements for the Torque particle J
//       dh_j_dQ_00[k] = u[0] *  arm_j[0][k];
//       dh_j_dQ_01[k] = u[1] *  arm_j[0][k];
//
//       dh_j_dQ_10[k] = u[0] *  arm_j[1][k];
//       dh_j_dQ_11[k] = u[1] *  arm_j[1][k];
//
//   }
//
//   // Helpers for exponential calculation
//   std::vector<double> G_hi_m1 = {G(h_i[0]-1.0), G(h_i[1]-1.0), G(h_i[2]-1.0)};
//   std::vector<double> G_hj_p1 = {G(h_j[0]+1.0), G(h_j[1]+1.0), G(h_j[2]+1.0)};
//
//   // optimiser helpers
//   double sum_G_hj_p1 = 0.0;
//   double G_r_rHB = G(r-r_HB);
//
//   double sum_dphiHB_dx = 0.0;
//   double sum_dphiHB_dy = 0.0;
//
//   double sum_dphiHB_dQ_j_00 = 0.0;
//   double sum_dphiHB_dQ_j_01 = 0.0;
//   double sum_dphiHB_dQ_j_10 = 0.0;
//   double sum_dphiHB_dQ_j_11 = 0.0;
//
//   // Assign the values from above
//   for(unsigned k=0;k<3;k++)
//     {
//         // calculate sum
//         sum_G_hj_p1 += G_hj_p1[k];
//
//       // h_i
//     //   h_i[k] =
// 	// ((*PartI->Q_pt(0,0)) * PartI->arm(0,k) + (*PartI->Q_pt(1,0)) * PartI->arm(1,k)) * u[0] +
// 	// ((*PartI->Q_pt(0,1)) * PartI->arm(0,k) + (*PartI->Q_pt(1,1)) * PartI->arm(1,k)) * u[1];
//       // h_i[k] = (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) * u[0]
//       //        + (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) * u[1];
//
//        // h_j
//      //   h_j[k] =
// 	 // ((*PartJ->Q_pt(0,0)) * PartJ->arm(0,k) + (*PartJ->Q_pt(1,0)) * PartJ->arm(1,k)) * u[0] +
// 	 // ((*PartJ->Q_pt(0,1)) * PartJ->arm(0,k) + (*PartJ->Q_pt(1,1)) * PartJ->arm(1,k)) * u[1];
//        // h_j[k] = (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) * u[0]
// 	   //        + (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) * u[1];
//
//      //   dG_dh_i[k] =
// 	 // -(h_i[k] - 1.0) / (Sigma_HB * Sigma_HB) * G(h_i[k] - 1.0);
//      //
//      //   // dG_dh_j
//      //   dG_dh_j[k] =
// 	 // -(h_j[k] + 1.0) / (Sigma_HB * Sigma_HB) * G(h_j[k] + 1.0);
//        dG_dh_i[k] = -(h_i[k] - 1.0) / (Sigma_HB * Sigma_HB) * G_hi_m1[k];
//
//        // dG_dh_j
//        dG_dh_j[k] = -(h_j[k] + 1.0) / (Sigma_HB * Sigma_HB) * G_hj_p1[k];
//
//
//      //   // dh_i_dx_i
//      // //   dh_i_dx_i[k] =
// 	 // // r_x * h_i[k] / (r * r)
// 	 // // - ((*PartI->Q_pt(0,0)) * PartI->arm(0,k) + (*PartI->Q_pt(1,0)) * PartI->arm(1,k)) / r;
//      //   dh_i_dx_i[k] = r_x * h_i[k] / (r * r)
//      //                - (Q_i[0][0] * arm_i[0][k] + Q_i[1][0] * arm_i[1][k]) / r;
//      //
//      //   // dh_j_dx_i
//      // //   dh_j_dx_i[k] =
// 	 // // r_x * h_j[k] / (r * r)
// 	 // // - ((*PartJ->Q_pt(0,0)) * PartJ->arm(0,k) + (*PartJ->Q_pt(1,0)) * PartJ->arm(1,k)) / r;
//      //    dh_j_dx_i[k] = r_x * h_j[k] / (r * r)
//      //                 - (Q_j[0][0] * arm_j[0][k] + Q_j[1][0] * arm_j[1][k]) / r;
//      //
//      //   // dh_i_dy_i
//      // //   dh_i_dy_i[k] =
// 	 // // r_y * h_i[k] / (r * r)
// 	 // // - ((*PartI->Q_pt(0,1)) * PartI->arm(0,k) + (*PartI->Q_pt(1,1)) * PartI->arm(1,k)) / r;
//      //    dh_i_dy_i[k] = r_y * h_i[k] / (r * r)
//      //                 - (Q_i[0][1] * arm_i[0][k] + Q_i[1][1] * arm_i[1][k]) / r;
//      //
//      //   // dh_j_dy_i
//      // //   dh_j_dy_i[k] =
// 	 // // r_y * h_j[k] / (r * r)
// 	 // // - ((*PartJ->Q_pt(0,1)) * PartJ->arm(0,k) + (*PartJ->Q_pt(1,1)) * PartJ->arm(1,k)) / r;
//      //    dh_j_dy_i[k] = r_y * h_j[k] / (r * r)
//      //                 - (Q_j[0][1] * arm_j[0][k] + Q_j[1][1] * arm_j[1][k]) / r;
//
//        // matrix elements for the Torque particle I
//        // dh_i_dQ_00[k] = u[0] *  PartI->arm(0,k);
//        // dh_i_dQ_01[k] = u[1] *  PartI->arm(0,k);
//        //
//        // dh_i_dQ_10[k] = u[0] *  PartI->arm(1,k);
//        // dh_i_dQ_11[k] = u[1] *  PartI->arm(1,k);
//        //
//        // // matrix elements for the Torque particle J
//        // dh_j_dQ_00[k] = u[0] *  PartJ->arm(0,k);
//        // dh_j_dQ_01[k] = u[1] *  PartJ->arm(0,k);
//        //
//        // dh_j_dQ_10[k] = u[0] *  PartJ->arm(1,k);
//        // dh_j_dQ_11[k] = u[1] *  PartJ->arm(1,k);
//        // dh_i_dQ_00[k] = u[0] *  arm_i[0][k];
//        // dh_i_dQ_01[k] = u[1] *  arm_i[0][k];
//        //
//        // dh_i_dQ_10[k] = u[0] *  arm_i[1][k];
//        // dh_i_dQ_11[k] = u[1] *  arm_i[1][k];
//        //
//        // // matrix elements for the Torque particle J
//        // dh_j_dQ_00[k] = u[0] *  arm_j[0][k];
//        // dh_j_dQ_01[k] = u[1] *  arm_j[0][k];
//        //
//        // dh_j_dQ_10[k] = u[0] *  arm_j[1][k];
//        // dh_j_dQ_11[k] = u[1] *  arm_j[1][k];
//
//        // calculate sum
//        sum_dphiHB_dx += dG_dh_j[k] * dh_j_dx_i[k];
//        sum_dphiHB_dy += dG_dh_j[k] * dh_j_dy_i[k];
//
//        sum_dphiHB_dQ_j_00 += dG_dh_j[k] * dh_j_dQ_00[k];
//        sum_dphiHB_dQ_j_01 += dG_dh_j[k] * dh_j_dQ_01[k];
//        sum_dphiHB_dQ_j_10 += dG_dh_j[k] * dh_j_dQ_10[k];
//        sum_dphiHB_dQ_j_11 += dG_dh_j[k] * dh_j_dQ_11[k];
//      }
//
//    // k over the arms
//   #pragma omp simd
//    for(unsigned k=0;k<3;k++)
//      {
//          // ---------------- HYDROGEN BOND LINEAR --------------- //
//          dphiHB_dx += sum_G_hj_p1 * (dG_dr * dr_dx * G_hi_m1[k]
//                                     + G_r_rHB * dG_dh_i[k] * dh_i_dx_i[k])
//                     + sum_dphiHB_dx * G_r_rHB * G_hi_m1[k];
//
//          dphiHB_dy += sum_G_hj_p1 * (dG_dr * dr_dy * G_hi_m1[k]
//                                     + G_r_rHB * dG_dh_i[k] * dh_i_dy_i[k])
//                     + sum_dphiHB_dy * G_r_rHB * G_hi_m1[k];
//
//         // --------------- HYDROGEN BOND TORQUE I -------------- //
//         dphiHB_dQ_i_00 += dG_dh_i[k] * dh_i_dQ_00[k] * sum_G_hj_p1;
//         dphiHB_dQ_i_01 += dG_dh_i[k] * dh_i_dQ_01[k] * sum_G_hj_p1;
//         dphiHB_dQ_i_10 += dG_dh_i[k] * dh_i_dQ_10[k] * sum_G_hj_p1;
//         dphiHB_dQ_i_11 += dG_dh_i[k] * dh_i_dQ_11[k] * sum_G_hj_p1;
//
//         // --------------- HYDROGEN BOND TORQUE J -------------- //
//         dphiHB_dQ_j_00 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_00;
//         dphiHB_dQ_j_01 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_01;
//         dphiHB_dQ_j_10 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_10;
//         dphiHB_dQ_j_11 +=  G_hi_m1[k] * sum_dphiHB_dQ_j_11;
//
//         // ---------------- HYDROGEN POTENTIAL --------------- //
//         *MolPt->potential_pt() += Epsilon_HB * G_r_rHB
//                                   * G_hi_m1[k] * sum_G_hj_p1;
//
//      //   for(unsigned l=0;l<3;l++)
// 	 // {
// 	 //   // ---------------- HYDROGEN BOND LINEAR --------------- //
// 	 //   // dphiHB_dx +=
// 	 //   //     Epsilon_HB * dG_dr     * dr_dx         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
// 	 //   //   + Epsilon_HB * G(r-r_HB) * dG_dh_i[k]    * dh_i_dx_i[k]  * G(h_j[l]+1.0)
// 	 //   //   + Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dx_i[l];
//      //   //
// 	 //   // dphiHB_dy +=
// 	 //   //     Epsilon_HB * dG_dr     * dr_dy         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
// 	 //   //   + Epsilon_HB * G(r-r_HB) * dG_dh_i[k]    * dh_i_dy_i[k]  * G(h_j[l]+1.0)
// 	 //   //   + Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dy_i[l];
//      //   // dphiHB_dx +=
//      //   //     dG_dr     * dr_dx         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
//      //   //   + G(r-r_HB) * dG_dh_i[k]    * dh_i_dx_i[k]  * G(h_j[l]+1.0)
//      //   //   + G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dx_i[l];
//      //   //
//      //   // dphiHB_dy +=
//      //   //      dG_dr     * dr_dy         * G(h_i[k]-1.0) * G(h_j[l]+1.0)
//      //   //   + G(r-r_HB) * dG_dh_i[k]    * dh_i_dy_i[k]  * G(h_j[l]+1.0)
//      //   //   +  G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l]    * dh_j_dy_i[l];
//      //
//      //     // dphiHB_dx += dG_dr * dr_dx * G_hi_m1[k] * G_hj_p1[l]
//      //     //   +  G_r_rHB * dG_dh_i[k]    * dh_i_dx_i[k]  * G_hj_p1[l]
//      //     //   +  G_r_rHB * G_hi_m1[k] * dG_dh_j[l]    * dh_j_dx_i[l];
//      //     //
//      //     // dphiHB_dy +=
//      //     //      dG_dr     * dr_dy         * G_hi_m1[k]  * G_hj_p1[l]
//      //     //   +  G_r_rHB * dG_dh_i[k]    * dh_i_dy_i[k]  * G_hj_p1[l]
//      //     //   +  G_r_rHB *G_hi_m1[k]  * dG_dh_j[l]    * dh_j_dy_i[l];
//      //
//      //   // dphiHB_dx += G_r_rHB * G_hi_m1[k] * dG_dh_j[l] * dh_j_dx_i[l];
//      //   //
//      //   // dphiHB_dy += G_r_rHB * G_hi_m1[k] * dG_dh_j[l] * dh_j_dy_i[l];
//      //
// 	 //   // --------------- HYDROGEN BOND TORQUE I -------------- //
// 	 //   // dphiHB_dQ_i_00 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * dG_dh_i[k] * dh_i_dQ_00[k] * G(h_j[l]+1.0);
//      //   //
// 	 //   // dphiHB_dQ_i_01 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * dG_dh_i[k] * dh_i_dQ_01[k] * G(h_j[l]+1.0);
//      //   //
// 	 //   // dphiHB_dQ_i_10 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * dG_dh_i[k] * dh_i_dQ_10[k] * G(h_j[l]+1.0);
//      //   //
// 	 //   // dphiHB_dQ_i_11 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * dG_dh_i[k] * dh_i_dQ_11[k] * G(h_j[l]+1.0);
//      //   // dphiHB_dQ_i_00 += dG_dh_i[k] * dh_i_dQ_00[k] * G(h_j[l]+1.0);
//      //   //
//      //   // dphiHB_dQ_i_01 += dG_dh_i[k] * dh_i_dQ_01[k] * G(h_j[l]+1.0);
//      //   //
//      //   // dphiHB_dQ_i_10 += dG_dh_i[k] * dh_i_dQ_10[k] * G(h_j[l]+1.0);
//      //   //
//      //   // dphiHB_dQ_i_11 += dG_dh_i[k] * dh_i_dQ_11[k] * G(h_j[l]+1.0);
//      //
// 	 //   // --------------- HYDROGEN BOND TORQUE J -------------- //
// 	 //   // dphiHB_dQ_j_00 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_00[l];
//      //   //
// 	 //   // dphiHB_dQ_j_01 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_01[l];
//      //   //
// 	 //   // dphiHB_dQ_j_10 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_10[l];
//      //   //
// 	 //   // dphiHB_dQ_j_11 +=
// 	 //   //   Epsilon_HB * G(r-r_HB) * G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_11[l];
//      //   // dphiHB_dQ_j_00 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_00[l];
//      //   //
// 	 //   // dphiHB_dQ_j_01 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_01[l];
//      //   //
// 	 //   // dphiHB_dQ_j_10 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_10[l];
//      //   //
// 	 //   // dphiHB_dQ_j_11 += G(h_i[k]-1.0) * dG_dh_j[l] * dh_j_dQ_11[l];
//      //
//      //   // dphiHB_dQ_j_00 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_00[l];
//      //   // dphiHB_dQ_j_01 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_01[l];
//      //   // dphiHB_dQ_j_10 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_10[l];
//      //   // dphiHB_dQ_j_11 +=  G_hi_m1[k] * dG_dh_j[l] * dh_j_dQ_11[l];
//      //
// 	 //   // // ---------------- HYDROGEN POTENTIAL --------------- //
// 	 //   // *MolPt->potential_pt() += Epsilon_HB * G(r - r_HB) * G(h_i[k]-1) * G(h_j[l]+1);
// 	 // }
//      }
//
// //    // add force for particle I
// // #pragma omp atomic
// //     *PartI->f_pt(0) -= dphiHB_dx;
// //
// // #pragma omp atomic
// //     *PartI->f_pt(1) -= dphiHB_dy;
// //
// //    // add force for particle J
// // #pragma omp atomic
// //     *PartJ->f_pt(0) += dphiHB_dx;
// //
// // #pragma omp atomic
// //     *PartJ->f_pt(1) += dphiHB_dy;
// //
// // // add up the forces between particles
// // forces[0] -= dphiHB_dx;
// // forces[1] -= dphiHB_dy;
//
//    // add force for particle I
// #pragma omp atomic
//     *PartI->f_pt(0) -= Epsilon_HB * dphiHB_dx;
//
// #pragma omp atomic
//     *PartI->f_pt(1) -= Epsilon_HB * dphiHB_dy;
//
//    // add force for particle J
// #pragma omp atomic
//     *PartJ->f_pt(0) += Epsilon_HB * dphiHB_dx;
//
// #pragma omp atomic
//     *PartJ->f_pt(1) += Epsilon_HB * dphiHB_dy;
//
// // add up the forces between particles
// forces[0] -= Epsilon_HB * dphiHB_dx;
// forces[1] -= Epsilon_HB * dphiHB_dy;
//
//
//    // add torque to overall torque
// // #pragma omp atomic
// //    *PartI->tau_pt(0) +=
// //      *PartI->Q_pt(0,1) * dphiHB_dQ_i_00 + *PartI->Q_pt(1,1) * dphiHB_dQ_i_10
// //      - *PartI->Q_pt(0,0) * dphiHB_dQ_i_01 - *PartI->Q_pt(1,0) * dphiHB_dQ_i_11;
// //
// // #pragma omp atomic
// //    *PartJ->tau_pt(0) +=
// //      *PartJ->Q_pt(0,1) * dphiHB_dQ_j_00 + *PartJ->Q_pt(1,1) * dphiHB_dQ_j_10
// //      - *PartJ->Q_pt(0,0) * dphiHB_dQ_j_01 - *PartJ->Q_pt(1,0) * dphiHB_dQ_j_11;
// #pragma omp atomic
//    *PartI->tau_pt(0) += Epsilon_HB * G_r_rHB * ( Q_i[0][1] * dphiHB_dQ_i_00
//                                                + Q_i[1][1] * dphiHB_dQ_i_10
//                                                - Q_i[0][0] * dphiHB_dQ_i_01
//                                                - Q_i[1][0] * dphiHB_dQ_i_11);
//
// #pragma omp atomic
//    *PartJ->tau_pt(0) += Epsilon_HB * G_r_rHB * ( Q_j[0][1] * dphiHB_dQ_j_00
//                                                + Q_j[1][1] * dphiHB_dQ_j_10
//                                                - Q_j[0][0] * dphiHB_dQ_j_01
//                                                - Q_j[1][0] * dphiHB_dQ_j_11);
//
//
//      // deleta arrays
//      // delete [] Q_i;
//      // delete [] Q_j;
//      // delete [] arm_i;
//      // delete [] arm_j;
//
//     // return the forces between the particles
//     return forces;
// };
//
//
// // function for computing the potential
// void MercBenz::compute_potential(Molecule* Mol_pt)
// {
//   printf("ERROR: This is not how the potential is calculated!");
//   exit(-1);
// };
//
//
// /******************************************************************************
//                           Lennard-Jones Cluster
//  *****************************************************************************/
// // Function which initialisese the system
// void LJcluster::system_initialise(Molecule* MolPt,
// 				  double& sigma,
// 				  double& epsilon,
// 				  const double& L_d)
// {
//   // check that the number of dimension is 2
//   if(2 != MolPt->dim())
//     {
//       printf("\nERROR: the dimension given to MercBenz is %i but should be 2\n\n",MolPt->dim());
//       exit(-1);
//     }
//
//   // set the sigma in the Lennard-Jones potential
//   Sigma_pt = &sigma;
//
//   // set the epsilon in the Lennard-Jones potential
//   Epsilon_pt = &epsilon;
//
//   unsigned number_of_particles = MolPt->nparticle();
//
//   // put particles into a perfect square
//   // for(unsigned i=0; i<number_of_particles; i++)
//   //   {
//   //       MolPt->particle_pt(i)->rigid_body()
//   //   }
// };
//
// // function for computing the force
// std::vector<double> LJcluster::compute_pair_force(Molecule* MolPt,
// 				   Particle* PartI,
// 				   Particle* PartJ,
// 				   const double& r,
// 				   const double& r_x,
// 				   const double& r_y)
// {
//   // tracker for pair force
//   std::vector<double> forces(2, 0.0);
//
//   // Lennard Jones Foxrce and Potential
//   forces = LJ_force(PartI,PartJ,r,r_x,r_y);
//
// #pragma omp atomic
//   *MolPt->potential_pt() += LJ_potential(r);
//
//   // return the pair force
//   return forces;
// };
//
//
// /******************************************************************************
//                             NN Polly System
//  *****************************************************************************/
// // Function for initialising the system
// void NNPolly::system_initialise(Molecule* MolPt,
// 				const double& L_d)
// {
//   // dereference the number of particles
//   NumPart = MolPt->nparticle();
//
//   double SepX = L_d/ ((double)NumPart+1.0);
//
//   // loop over all the particles and put them into the grid
//   for(unsigned i=0;i<NumPart;i++)
//     {
//       *MolPt->particle_pt(i)->q_pt(0) = 0.5 * SepX + i * SepX;
//     }
// };
//
// // function for computing the potential
// void NNPolly::compute_potential(Molecule* Mol_pt)
// {
// };
//
//
// /******************************************************************************
//                     Hydrogen Spin system
//  *****************************************************************************/
// // initialise the system
// void HSpin::system_initialise(Molecule* molecule_pt,
//                               const double& e_field)
// {
//     // set the strength of the field
//     external_field = e_field;
//
//     // compute the force
//     this->compute_force(molecule_pt);
// }
//
// // function for computing the potential
// void HSpin::compute_potential(Molecule* Mol_pt)
// {
//   //empty
// }
//
// // compute force
// void HSpin::compute_force(Molecule* molecule_pt)
// {
//     // dereference helpers
//     double* f = NULL;
//     double* theta = NULL;
//
//     unsigned K = molecule_pt->nparticle();
//
//     // reset the magnitisation
//     magnetisation = 0.0;
//
//     // compute the magnitisation
//     for(unsigned i=0; i<K; i++)
//     {
//         // get the rotation of the particle
//         theta = molecule_pt->particle_pt(i)->q_pt(0);
//
//         // mod the particle rotation
//         *theta = fmod(*theta, 2.0*M_PI);
//
//         // calculate the magnitisation
//         magnetisation += cos(*theta);
//     }
//
//     // devide by numparts
//     magnetisation /= (double)K;
//
//     // assign the force on each particle
//     for(unsigned i=0; i<K; i++)
//     {
//         // dereference the particle
//         theta = molecule_pt->particle_pt(i)->q_pt(0);
//         f = molecule_pt->particle_pt(i)->f_pt(0);
//
//         // force is negative of the gradient
//         *f = -( magnetisation + external_field ) * sin(*theta);
//     }
//
//     // calculate the potential too
//     *molecule_pt->potential_pt() = -(double)K * ( 0.5 * pow(magnetisation, 2.0)
//                                            + external_field * magnetisation);
// };
