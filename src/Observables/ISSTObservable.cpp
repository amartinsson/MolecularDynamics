#include "ISSTObservable.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
IsstObservable::IsstObservable(InfiniteSwitchSimulatedTempering* method,
    const int& recf, const int& rect)
        : SystemObservable(recf, rect), method(method)
{
    // empty
}
