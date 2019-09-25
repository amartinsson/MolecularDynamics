#include "DISObservable.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
DisObservable::DisObservable(DoubleInfiniteSwitch* method,
    const int& recf, const int& rect)
        : SystemObservable(recf, rect), method(method)
{
    // empty
}
