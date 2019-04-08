#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include<vector>
#include<iostream>
// #include<omp.h>
// #include<stdlib.h>
// #include<math.h>
// #include<gsl/gsl_integration.h>
// #include<stdio.h>
// #include<gsl/gsl_rng.h>
// #include<gsl/gsl_randist.h>

#include "Molecules.hpp"
#include "Array.hpp"

using namespace::std;

class System
{
public:
	System(){};
	// calculate the force
	virtual void compute_force(Molecule* molecule_pt) = 0;
	// calculate the pair force
	virtual Vector compute_force(Molecule* molecule_pt,
										   		   Particle* particle_i,
				                           		   Particle* particle_j,
										   		   const double& r,
												   const Vector& dr) = 0;
};

// class Systems
// {
// protected:
// // Lennard Jones parameter sigma
// double* Sigma_pt;
//
// // Lennard Jones parameter epsilon
// double* Epsilon_pt;
//
// public:
//
// //constructor
// Systems();
//
// // function for the Lennard Jones potential
// double LJ_potential(const double& r);
//
// // function for the derivative of Lennard Jones potential
// std::vector<double> LJ_force(Particle* PartI,
// 		Particle* PartJ,
// 		const double& r,
// 		const double& r_x,
// 		const double& r_y);
//
//   // function for initilasing the system
//   void system_initialise(std::vector<Molecule*>& Mol_pt,
// 			 const unsigned& nCopies);
//
//   // function for computing the force
//   void compute_force(Molecule* Mol_pt);
//
//   // function for computing the pair force
//   virtual std::vector<double> compute_pair_force(Molecule* MolPt,
// 				  Particle* PartI,
// 				  Particle* PartJ,
// 				  const double& r,
// 				  const double& r_x,
// 				  const double& r_y) = 0;
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
//
//   // function for computing the normalising constant
//   //  virtual double get_normalising_constant(const double& beta) = 0;
//
//   // function for determining region
//   void determine_region(std::vector<Molecule*>& Mol_pt);
//
//   // function which returns the average beta in ISST
//   double get_beta_bar(const double& V,
// 		      const unsigned& D,
// 		      const double& betaOne,
// 		      const double& betaTwo);
//
//   // function determining partition function Z(beta) if it exists
//   double domega_dbeta_over_omega(const double& beta, const unsigned& D);
//
//   virtual double get_exact_pdf(const double& x,
// 			       const double& beta){return 0.0;};
//
//   virtual double get_exact_q2(const double& beta) {return 0.0;};
//
//   virtual double get_exact_gauss(const double& beta,
// 				 const double& mu,
// 				 const double& sigma) {return 0.0;};
//
//   virtual std::vector<double> exact_pdf()
//   {
//     std::vector<double> bajs;
//     return bajs;
//   }
//
// };

// /******************************************************************************
//                            An-Harmonic Oscillator
//  *****************************************************************************/
// #include <fstream>
//
// class AnHarmonicOscillator : public Systems
// {
// private:
//
//   // double holding the epsilon parameter
//   static double Epsilon;
//   static double Eta;
//
//   static double gauss_mu;
//   static double gauss_sigma;
//
// public:
//
//   // constructor - empty
//   AnHarmonicOscillator(){};
//
//   // initialise the system
//   void system_initialise
//   (Molecule* MolPt,
//    const double& epsilon,
//    const double& eta);
//
//   // compute the force using the molecule
//   void compute_force
//   (Molecule* MolPt);
//
//   // compute the force using particle
//   void compute_force
//   (Particle* PartPt,
//    const unsigned& DIM);
//
//   // compute the potential
//   void compute_potential
//   (Molecule* MolPt);
//
//   // function for computing the pair force
// 	// function for computing the pair force
// 	std::vector<double> compute_pair_force(Molecule* MolPt,
// 					Particle* PartI,
// 					Particle* PartJ,
// 					const double& r,
// 					const double& r_x,
// 					const double& r_y) {return std::vector<double>(1,0.0);}
//
//   // function returning the partition function integrand at position x
//   static double z_beta(double x, void * params);
//
//   // function returning the partition function integrand at position x
//   static double q2_beta(double x, void * params);
//
//   // function returning the partition function integrand at position x
//   static double gauss_beta(double x, void * params);
//
//   // function which returns the exact pdf - calculated numerically
//   double get_exact_pdf(const double& x,
// 		       const double& Beta);
//
//   // get the exact solution to < q^2 >
//   double get_exact_q2(const double& beta);
//
//   // get the exact solution to < q^2 >
//   double get_exact_gauss(const double& beta,
// 			 const double& mu,
// 			 const double& sigma);
// };
//
// /******************************************************************************
//                                Double Well
//  *****************************************************************************/
//
// class DoubleWell : public Systems
// {
//   // Double well system can be found at:
//   // http://www.ergodic.org.uk/MDdotM/MDdotM.php?Model=Double_Well&headerbar=1
//
//  private:
//
//   // Parameters as described on the website
//   static double a;
//   static double b;
//
//  public:
//
//   // Constructor
//   DoubleWell();
//
//   // function for initilasing the system
//   void system_initialise(Molecule* Mol_pt,
// 			 double& a,
// 			 double& b);
//
//   // function for computing the force
//   void compute_force(Molecule* Mol_pt);
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
//
//   // function returning the partition function integrand at position x
//   static double z_beta(double x, void * params);
//
//   // returns the normalisation constant at temperature beta
//   double get_normalisation_constant(const double& beta);
//
//   // return exact pdf
//   double get_exact_pdf(const double& x,
// 		       const double& beta);
//
// 					 // function for computing the pair force
// 				 	std::vector<double> compute_pair_force(Molecule* MolPt,
// 				 				  Particle* PartI,
// 				 				  Particle* PartJ,
// 				 				  const double& r,
// 				 				  const double& r_x,
// 				 				  const double& r_y) {return std::vector<double>(1,0.0);}
// };
//
// /******************************************************************************
//                                Trimer
//  *****************************************************************************/
//
// class Trimer : public Systems
// {
//   // Trimer system can be found at:
//   // http://www.ergodic.org.uk/MDdotM/MDdotM.php?Model=Trimer&headerbar=1
//
//  public:
//
//   // Constructor
//   Trimer();
//
//   // function for initilasing the system
//   void system_initialise(Molecule* Mol_pt,
// 			 const double& sigma,
// 			 const double& epsilon);
//
//   // function for computing the force
//   void compute_force(Molecule* Mol_pt);
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
// };
//
// /******************************************************************************
//                                   Vanden
//  *****************************************************************************/
// class Vanden : public Systems
// {
//   // Vanden System is an implementation of a system as described in:
//   // (18) from DOI: 10.1063/1.4790706
//
//  private:
//
//   // holder for Lambda constant
//   double Lambda;
//
//   // R_k vector as defined in the Inifinte Replica Exchange paper by Jianfeng
//   // and Eric
//   double R;
//
//  public:
//
//   // Constructor
//   Vanden();
//
//   // Function for initilasing the system
//   void system_initialise(Molecule* Mol_pt,
// 			 double& lambda);
//
//   // Function that return the force scaling in Infinite Replica Exchange R_k
//   double* r_pt();
//
//   // Function for computing the force
//   void compute_force(Molecule* Mol_pt);
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
//
//   // function for computing the normalising constant
//   double get_normalising_constant(const double& beta){return -1.0;};
// };
//
// /******************************************************************************
//                                 Mercedes-Benz
//  *****************************************************************************/
// class MercBenz : public Systems
// {
//   // Mercedes Benz system is an implemntation of a system described in:
//   // https://arxiv.org/pdf/1304.3232.pdf
//
//  private:
//
//   // Total Number of Particles
//   unsigned NumPartRow;
//   unsigned NumPartCol;
//
//   // The Total Box is [-Lx, Lx, -Ly, Ly]
//   double Lx;
//   double Ly;
//
//   unsigned test;
//
//   // Holder for the cutoff in the force calculation
//   double CutOff;
//
//   // Holder for the number of boxes
//   unsigned TotNumBoxes;
//
//  // sigma in the Hydrogen Bond part
//   double Sigma_HB;
//
//   // radius Hydrogen Bond part
//   double r_HB;
//
//   // epsilon Hydrogen Bond part
//   double Epsilon_HB;
//
//   // Gaussian function
//   double G(const double& x);
//
//  public:
//
//   // empty constructor
//   MercBenz();
//
//   // Function for initialising the system
//   void system_initialise(Molecule* MolPt,
// 			 double& sigma,
// 			 double& epsilon,
// 			 double& sigma_hb,
// 			 double& r_hb,
// 			 double& epsilon_hb,
// 			 const double& L_d);
//
//   // function for computing the force
//   std::vector<double> compute_pair_force(Molecule* MolPt,
// 			  Particle* PartI,
// 			  Particle* PartJ,
// 			  const double& r,
// 			  const double& r_x,
// 			  const double& r_y);
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
// };
//
// /******************************************************************************
//                           Lennard-Jones Cluster
//  *****************************************************************************/
// class LJcluster : public Systems
// {
//  public:
//
//   // empty constructor
//   LJcluster(){};
//
//   // Function for initialising the system
//   void system_initialise(Molecule* MolPt,
// 			 double& sigma,
// 			 double& epsilon,
// 			 const double& L_d);
//
//   // function for computing the pair force
// 	// function for computing the pair force
// 	std::vector<double> compute_pair_force(Molecule* MolPt,
// 					Particle* PartI,
// 					Particle* PartJ,
// 					const double& r,
// 					const double& r_x,
// 					const double& r_y);
// };
//
// /******************************************************************************
//                     Nearest Neighbour Polymere
//  *****************************************************************************/
// class NNPolly : public Systems
// {
//    private:
//
//   // Total Number of Particles
//   unsigned NumPart;
//
//   // box length
//   double L;
//
//  public:
//
//   // empty constructor
//   NNPolly(){};
//
//   // Function for initialising the system
//   void system_initialise(Molecule* MolPt,
// 			 const double& L_d);
//
//   // function for computing the force
// 	// function for computing the pair force
// 	std::vector<double> compute_pair_force(Molecule* MolPt,
// 					Particle* PartI,
// 					Particle* PartJ,
// 					const double& r,
// 					const double& r_x,
// 					const double& r_y){return std::vector<double>(1,0.0);}
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
// };
//
//
// /******************************************************************************
//                     Heisenburg Spin system
//  *****************************************************************************/
// class HSpin : public Systems
// {
//    private:
//
//   // total hold the magnitisation
//   double magnetisation;
//
//   // holder for the external field
//   double external_field;
//
//  public:
//   // empty constructor
//   HSpin(){};
//
//   // initialise the system
//   void system_initialise(Molecule* MolPt, const double& H);
//
//   // function for computing the potential
//   void compute_potential(Molecule* Mol_pt);
//
//   // compute force
//   void compute_force(Molecule* MolPt);
//
//   // compute pair force
// 	// function for computing the pair force
// 	std::vector<double> compute_pair_force(Molecule* MolPt,
// 					Particle* PartI,
// 					Particle* PartJ,
// 					const double& r,
// 					const double& r_x,
// 					const double& r_y){return std::vector<double>(1,0.0);}
//
//   // return the magnitisation
//   double magnitisation()
//   {
//     return magnetisation;
//   }
//
// };


#endif
