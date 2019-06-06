#include "SystemVolume.hpp"

using namespace::std;

// constructor
SystemVolume::SystemVolume(const Matrix& s, const unsigned& n,
                               const int& recf, const int& rect)
        : SystemObservable(recf, rect), S(&s), N(n)
{
    Volume = new AverageObservable();
}

// destructor
SystemVolume::~SystemVolume()
{
    delete Volume;
}

// instant
double SystemVolume::get_instant()
{
    return Volume->get_instant();
}

// average
double SystemVolume::get_average()
{
    return Volume->get_average();
}

// update function
void SystemVolume::update()
{
    if(recStep())
    {
        // make system observation
        Volume->observe((*S).det());
    }
}
