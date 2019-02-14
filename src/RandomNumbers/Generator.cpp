#include "Generator.hpp"

using namespace std;

/******************************************************************************
                    Normal Random Number Generator Class
 *****************************************************************************/
double NormalGenerator::operator ()()
{
    return this->distribution(Generator::generator);
}

/******************************************************************************
                    Uniform Random Number Generator Class
 *****************************************************************************/
double UniformGenerator::operator ()()
{
    if(&real_distribution != NULL)
        return this->real_distribution(Generator::generator);
    else
        return this->int_distribution(Generator::generator);
}
