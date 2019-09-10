#include "ReplicaCopy.hpp"

using namespace std;

ReplicaCopy::ReplicaCopy(Molecule* original, const unsigned& Ncopies)
    : Ncopies(Ncopies)
{

    // Singelton(const Vector& q_0, const Vector& p_0, const Matrix& m,
    //           const double& kt, const unsigned& dim);
    Vector q0 = original->particle(0).q;
    Vector p0 = original->particle(0).p;
    Matrix m  = original->particle(0).m;
    double kt = original->kt();
    double dim = original->dim();

    for(unsigned i=0; i<Ncopies; i++) {
        if(i == 0) {
            molecule.push_back(original);
        }
        else {
            molecule.push_back(new Singelton(q0, p0, m, kt, dim));
        }
    }
}

ReplicaCopy::~ReplicaCopy()
{
    for(unsigned i=0; i<Ncopies; i++) {
        delete molecule[i];
    }
}

Molecule* ReplicaCopy::operator() (const unsigned& i) const
{
    return molecule[i];
}

unsigned ReplicaCopy::size() const
{
    return Ncopies;
}

void ReplicaCopy::exchange_pos_mom_tem(Molecule* m1, Molecule* m2)
{
    // position
    Vector q1 = m1->particle(0).q;
    Vector q2 = m2->particle(0).q;

    m1->particle(0).q = q2;
    m2->particle(0).q = q1;

    // momentum
    Vector p1 = m1->particle(0).p;
    Vector p2 = m2->particle(0).p;

    m1->particle(0).p = p2;
    m2->particle(0).p = p1;

    // temperature
    double b1 = m1->beta();
    double b2 = m2->beta();

    m1->set_beta(b2);
    m2->set_beta(b1);
}

void ReplicaCopy::exchange_pos_tem(Molecule* m1, Molecule* m2)
{
    // position
    Vector q1 = m1->particle(0).q;
    Vector q2 = m2->particle(0).q;

    m1->particle(0).q = q2;
    m2->particle(0).q = q1;

    // temperature
    double b1 = m1->beta();
    double b2 = m2->beta();

    m1->set_beta(b2);
    m2->set_beta(b1);
}
