#include "SystemPressure.hpp"

using namespace::std;

// constructor
SystemPressure::SystemPressure(const Matrix& s, const Matrix& v,
                               const Matrix& k,
                               const int& recf, const int& rect)
        : SystemObservable(recf, rect), S(&s), V(&v), K(&k)
{
    Pressure = new AverageObservable();
}

// destructor
SystemPressure::~SystemPressure()
{
    delete Pressure;
}

// instant
double SystemPressure::get_instant()
{
    return Pressure->get_instant();
}

// average
double SystemPressure::get_average()
{
    return Pressure->get_average();
}

// update function
void SystemPressure::update()
{
    if(recStep())
    {
        double p = calculate_pressure(*S, *V, *K);

        // make system observation
        Pressure->observe(p);
    }
}

// // print function at timeindex
// void SystemPressure::print(const char* file_name,
//                            const double& time,
//                            const unsigned& index)
// {
//     // convert the filename to a string
//     std::string name(file_name);
//
//     // open the file to write to
//     char filename[50];
//     sprintf(filename, "Observables/%s_%d.csv", (name).c_str(), index);
//     FILE* file = fopen(filename, "a");
//
//     fprintf(file, "%1.7e, %1.7e, %1.7e\n",
//             time, this->get_instant(), this->get_average());
//
//     // close the file
//     fclose(file);
// }
//
// // print function at time
// void SystemPressure::print(const char* file_name,
//                            const double& time)
// {
//     // convert the filename to a string
//     std::string name(file_name);
//
//     // open the file to write to
//     char filename[50];
//     sprintf(filename, "Observables/%s.csv", (name).c_str());
//     FILE* file = fopen(filename, "a");
//
//     fprintf(file, "%1.7e, %1.7e, %1.7e\n",
//             time, this->get_instant(), this->get_average());
//
//     // close the file
//     fclose(file);
// }

// function which calculates the pressure
double SystemPressure::calculate_pressure(const Matrix& s,
                                          const Matrix& v,
                                          const Matrix& k)
{
    // pressure
    double pressure = 0.0;

    if(s.size()[0] > 2)
        pressure = -1.0 / (3.0 * s(1,1) * s(2,2)) * (k(0,0) + v(0,0))
                   -1.0 / (3.0 * s(0,0) * s(2,2)) * (k(1,1) + v(1,1))
                   -1.0 / (3.0 * s(0,0) * s(1,1)) * (k(2,2) + v(2,2));
    else
        pressure = -1.0 / (2.0 * s(1,1)) * (k(0,0) + v(0,0))
                   -1.0 / (2.0 * s(0,0)) * (k(1,1) + v(1,1));

    return pressure;
}
