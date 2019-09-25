#include "ISObservable.hpp"

using namespace::std;

/******************************************************************************
                Abstract SimulatedTempering Observable Base Class
 *****************************************************************************/
IsObservable::IsObservable(InfiniteSwitch* method,
    const int& recf, const int& rect)
        : SystemObservable(recf, rect), method(method)
{
    // empty
}
