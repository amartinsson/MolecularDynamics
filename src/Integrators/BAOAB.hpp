#include "Langevin.hpp"
#include "System.hpp"

using namespace::std;

/******************************************************************************
                                 BAOAB Class
 *****************************************************************************/
class BAOAB : public Langevin
{
public:
    BAOAB(const double& beta, const double& gamma, const double& gamma_rot,
          const double& time_step, System* system_pt, const int& seed);
    // destructor
    ~BAOAB();
    // must have integrate
    void integrate(Molecule* molecule_pt);
    // set with npt grid
    void integrate_with_npt_grid(const Matrix& Szero, const double& cut_off,
                                 Molecule* molecule_pt, const double& mass,
                                 const double& target_press,
                                 const double& gamma_box, const int& recf,
                                 const int& rect);
    // set which npt integration version
    void set_npt_integrator_version(const unsigned& version);

private:
    double Time_Step;
    unsigned Step;
    unsigned Npt_version;

    bool With_npt;
    bool With_st;

    void nvt_integration(Molecule* molecule_pt);
    void npt_integration(Molecule* molecule_pt);
};
