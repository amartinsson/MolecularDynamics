#include "Generator.hpp"

using namespace std;

/******************************************************************************
                    Normal Random Number Generator Class
 *****************************************************************************/
double NormalGenerator::operator() ()
{
    return this->distribution(Generator::generator);
}

/******************************************************************************
                    Uniform Random Number Generator Class
 *****************************************************************************/
double UniformGenerator::operator() ()
{
    return this->real_distribution(Generator::generator);
}

/******************************************************************************
                Exponential Random Number Generator Class
 *****************************************************************************/
double ExponentialGenerator::operator() ()
{
    return this->distribution(Generator::generator);
}

// return exponential with generator lambda
double ExponentialGenerator::operator() (const double& lambda) {
    return -log(uniform()) / lambda;
}
