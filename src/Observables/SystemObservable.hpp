#ifndef SYSTEMOBSERVABLE_HPP
#define SYSTEMOBSERVABLE_HPP

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class SystemObservable
{
public:
    // empty constructor
    SystemObservable(){};
    // empty destructor
    ~SystemObservable(){};
    // virtual function for updating
    virtual void update() = 0;
    // print function at time_index
    virtual void print(const char* file_name, const double& time,
                       const unsigned& index) = 0;
    // print function at time
    virtual void print(const char* file_name, const double& time) = 0;
};



#endif
