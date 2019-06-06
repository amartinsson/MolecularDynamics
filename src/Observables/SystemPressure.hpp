#ifndef SYSTEMPRESSURE_HPP
#define SYSTEMPRESSURE_HPP

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "SystemObservable.hpp"
#include "Molecules.hpp"

using namespace::std;

/******************************************************************************
                        System Temperature Class
 *****************************************************************************/
class SystemPressure : public SystemObservable
{
public:
    // constructor
    SystemPressure(const Matrix& s, const Matrix& v, const Matrix& k,
                   const int& recf, const int& rect);
    // destructor
    ~SystemPressure();
    // instant
    double get_instant();
    // average
    double get_average();
    // update function
    void update();

private:
    // private momentum temperature
    AverageObservable* Pressure;
    // Matrix pointers
    const Matrix* S;
    const Matrix* V;
    const Matrix* K;

    // calulate the pressure
    double calculate_pressure(const Matrix& s, const Matrix& v,
                              const Matrix& k);
};

#endif
