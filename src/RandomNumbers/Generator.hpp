#ifndef GENERATOR_HPP
#define GENERATOR_HPP
#include <random>

using namespace std;

/******************************************************************************
                    Base Random Number Generator Class
 *****************************************************************************/
class Generator
{
public:
    Generator(const int& seed) : generator(seed) {};
    virtual double operator ()() = 0;
    ~Generator() {};

protected:
    std::mt19937_64 generator;
};

/******************************************************************************
                    Normal Random Number Generator Class
 *****************************************************************************/
class NormalGenerator : public Generator
{
public:
    NormalGenerator(const double& mean, const double& std, const int& seed):
    Generator(seed), distribution(mean, std)
    {};

    // destructor
    ~NormalGenerator(){};

    // return normal number
    double operator ()();

private:
    std::normal_distribution<double> distribution;
};

/******************************************************************************
                    Uniform Random Number Generator Class
 *****************************************************************************/
class UniformGenerator : public Generator
{
public:
    UniformGenerator(const double& min, const double& max, const int& seed):
    Generator(seed), real_distribution(min, max) {};

    // destructor
    ~UniformGenerator(){};

    // return uniform real number
    double operator ()();

private:
    std::uniform_real_distribution<double> real_distribution;
};

#endif
