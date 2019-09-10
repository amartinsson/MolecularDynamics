#ifndef REPLICACOPY_HPP
#define REPLICACOPY_HPP

#include <Molecules.hpp>

using namespace std;

class ReplicaCopy
{
public:
    ReplicaCopy(Molecule* original, const unsigned& Ncopies);

    ~ReplicaCopy();

    void exchange_pos_mom_tem(Molecule* m1, Molecule* m2);
    void exchange_pos_tem(Molecule* m1, Molecule* m2);
    unsigned size() const;

    Molecule* operator() (const unsigned& i) const;

private:
    unsigned Ncopies;
    vector<Molecule*> molecule;
};
#endif
