#ifndef SYSTEMTRAJECTORY_HPP
#define SYSTEMTRAJECTORY_HPP

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
private:
    // Pointer to molecule object
    const Molecule* molecule_pt;
    // constant character object
    const char* file_name;
    // Pointer to simulation box object
    const Matrix* S;
};

#endif
