#ifndef SYSTEMOBSERVABLE_HPP
#define SYSTEMOBSERVABLE_HPP

#include <fstream>
#include <sstream>

using namespace::std;

/******************************************************************************
                        Abstract System Observable Base Class
 *****************************************************************************/
class SystemObservable
{
public:
    // empty constructor
    SystemObservable(const int& recf, const int& rect)
        : RecFreq(recf), RecCount(0), RecTotCount(0), RecThresh(rect) {};
    // empty destructor
    ~SystemObservable(){};
    // virtual function for updating
    virtual void update() = 0;
    // print function at time_index
    void print(const char* file_name, const double& time,
                       const unsigned& index);
    // print function at time
    void print(const char* file_name, const double& time);
    // get the average and instant vaalues
    virtual double get_average()=0;
    // get the instant value
    virtual double get_instant()=0;
    // reset
    void reset(const int& Recf, const int& Rect)
    {
        RecFreq = Recf;
        RecThresh = Rect;

        RecCount = 0;
        RecTotCount = 0;
    }


protected:
    // check that step is a record step
    bool recStep();

private:
    // recording frequency
    int RecFreq;
    int RecCount;

    // recording threshold
    int RecThresh;
    int RecTotCount;
};



#endif
