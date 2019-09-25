#ifndef SYSTEMTRAJECTORY_HPP
#define SYSTEMTRAJECTORY_HPP

#include "HistObservable.hpp"
#include "Molecules.hpp"

using namespace::std;

/******************************************************************************
                        System Trajectory Class
 *****************************************************************************/
class SystemTrajectory
{
public:
    // constructor
    SystemTrajectory(const Molecule* molecule_pt);
    // destructor
    ~SystemTrajectory() {/* empty */};
    // set the simbox to record Trajectory
    void set_simbox(const Matrix& s);
    // print function at time_index
    void print_positions(const char* file_name, const unsigned& index);
    // print simbox
    void print_simbox(const char* file_name, const unsigned& index);
    // print function at time_index
    void append_positions(const char* file_name, const double& time);
    void append_positions(const char* file_name, const unsigned& index,
        const double& time);
protected:
    // Pointer to molecule object
    const Molecule* molecule_pt;
    // constant character object
    const char* file_name;
    // Pointer to simulation box object
    const Matrix* S;
};

class SystemHistogramTrajectory : public SystemTrajectory
{
public:
    // constructor
    SystemHistogramTrajectory(const Molecule* molecule_pt,
        const vector<double>& min, const vector<double>& max,
            const vector<int>& N);
    // destructor
    ~SystemHistogramTrajectory();
    // update the hisgoram
    void update();
    void update(const double& weight);
    // print the final histogram
    void print(const char* file_name);
    void print(const char* file_name, const unsigned& index);

private:
    vector<HistObservable*> posDist;
};

#endif
