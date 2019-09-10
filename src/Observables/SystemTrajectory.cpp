#include "SystemTrajectory.hpp"

using namespace::std;

// constructor
SystemTrajectory::SystemTrajectory(const Molecule* molecule_pt)
            : molecule_pt(molecule_pt) {/* empty */}

// set the box object
void SystemTrajectory::set_simbox(const Matrix& s) { S = &s;}

// print function at time_index
void SystemTrajectory::print_positions(const char* file_name,
    const unsigned& index) {
        // open the file to write to
        char filename[50];
        sprintf(filename, "Observables/Frames/%s_%i.csv", file_name,
                index);
        FILE* file = fopen(filename, "w");

        for(auto& particle : molecule_pt->Particles) {
            // print all translational degrees of freedom
            for(unsigned i=0; i<particle.second->q.size(); i++) {
                if(i == particle.second->q.size()-1)
                    fprintf(file, "%.4f", particle.second->q(i));
                else
                    fprintf(file, "%.4f, ", particle.second->q(i));
            }

            // print the rotational degrees of freedom
            if(particle.second->rigid_body())
                fprintf(file, ", %.10f, %.10f, %.10f, %.10f\n",
                particle.second->Q(0,0), particle.second->Q(0,1),
                particle.second->Q(1,0), particle.second->Q(1,1));
            else
                fprintf(file, "\n");
        }

        // close the file
        fclose(file);
}

// print the simbox trajectory
void SystemTrajectory::print_simbox(const char* file_name,
    const unsigned& index) {
        // open the file to write to
        char filename[50];
        sprintf(filename, "Observables/Frames/%s_%i.csv", file_name,
                index);
        FILE* file = fopen(filename, "w");

        // loop over to print the simbox
        for(int i=0;i<10;i++)
        {
            if(i==0 | i==4)
              fprintf(file, "%.3f, %.3f, %.3f\n", 0.0, 0.0, 0.0);
            else if(i==1)
              fprintf(file, "%.3f, %.3f, %.3f\n", ((*S))(0,0), 0.0, 0.0);
            else if(i==2)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,0)+(*S)(0,1),
              (*S)(1,1), 0.0);
            else if(i==3)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,1), (*S)(1,1), 0.0);
            else if(i==5 || i==9)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,2), (*S)(1,2),
              (*S)(2,2));
            else if(i==6)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,0)+(*S)(0,2),
              (*S)(1,2), (*S)(2,2));
            else if(i==7)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,0)+(*S)(0,1)+(*S)(0,2),
              (*S)(1,1)+(*S)(1,2),(*S)(2,2));
            else if(i==8)
              fprintf(file, "%.3f, %.3f, %.3f\n", (*S)(0,1)+(*S)(0,2),
              (*S)(1,1)+(*S)(1,2), (*S)(2,2));
        }

        fclose(file);
    }

// append function at time_index
void SystemTrajectory::append_positions(const char* file_name,
    const double& time) {
        // open the file to write to
        char filename[50];
        sprintf(filename, "Observables/%s.csv", file_name);
        FILE* file = fopen(filename, "a");

        fprintf(file, "%1.7e", time);

        for(auto& particle : molecule_pt->Particles) {
            // print all translational degrees of freedom
            for(unsigned i=0; i<particle.second->q.size(); i++)
                fprintf(file, ", %.4f", particle.second->q(i));

            // print the rotational degrees of freedom
            if(particle.second->rigid_body())
                fprintf(file, ", %.10f, %.10f, %.10f, %.10f\n",
                particle.second->Q(0,0), particle.second->Q(0,1),
                particle.second->Q(1,0), particle.second->Q(1,1));
        }

        // end with new line

        fprintf(file, "\n");
        // close the file
        fclose(file);
    }

// append function at time_index
void SystemTrajectory::append_positions(const char* file_name,
    const unsigned& index, const double& time) {
        // open the file to write to
        char filename[50];
        sprintf(filename, "Observables/%s_%i.csv", file_name, index);
        FILE* file = fopen(filename, "a");

        fprintf(file, "%1.7e", time);

        for(auto& particle : molecule_pt->Particles) {
            // print all translational degrees of freedom
            for(unsigned i=0; i<particle.second->q.size(); i++)
                fprintf(file, ", %.4f", particle.second->q(i));

            // print the rotational degrees of freedom
            if(particle.second->rigid_body())
                fprintf(file, ", %.10f, %.10f, %.10f, %.10f\n",
                particle.second->Q(0,0), particle.second->Q(0,1),
                particle.second->Q(1,0), particle.second->Q(1,1));
        }

        // end with new line

        fprintf(file, "\n");
        // close the file
        fclose(file);
    }
