#ifndef SYSTEMVOLUME_HPP
#define SYSTEMVOLUME_HPP

#include "AverageObservable.hpp"
#include "Array.hpp"
#include "SystemObservable.hpp"

using namespace::std;

/******************************************************************************
                        System Temperature Class
 *****************************************************************************/
class SystemVolume : public SystemObservable
{
public:
    // constructor
    SystemVolume(const Matrix& s, const unsigned& n, const int& recf,
                const int& rect);
    // destructor
    ~SystemVolume();
    // instant
    double get_instant();
    // average
    double get_average();
    // update function
    void update();
    // reset
    void reset(const int& recf, const int& rect) {
        SystemObservable::reset(recf, rect);

        Volume->clear();
    }

private:
    // private momentum temperature
    AverageObservable* Volume;
    // Matrix pointers
    const Matrix* S;
    const unsigned N;
};

#endif
