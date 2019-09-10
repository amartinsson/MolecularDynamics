#include "ReplicaExchangeObservable.hpp"

using namespace::std;

/******************************************************************************
                Trajectory Replica Exchange Observable
 *****************************************************************************/
RETrajectoryObservable::RETrajectoryObservable(ReplicaExchange* method,
    const int& recf, const int& rect)
        : N(method->replica.size()), SystemObservable(recf, rect)
{
    for(unsigned i=0; i<N; i++) {
        trajectory.push_back(new SystemTrajectory(method->replica(i)));
    }
}

void RETrajectoryObservable::update()
{
    // do nothing
}


void RETrajectoryObservable::append_positions(const char* file_name,
    const double& time)
{
    for(unsigned i=0; i<N; i++) {
        trajectory[i]->append_positions(file_name, i, time);
    }
}
