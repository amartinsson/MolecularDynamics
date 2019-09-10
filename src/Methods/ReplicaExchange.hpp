#ifndef REPLICAEXCHANGE_HPP
#define REPLICAEXCHANGE_HPP

#include "Molecules.hpp"
#include "ReplicaCopy.hpp"
#include "Generator.hpp"

using namespace std;

class ReplicaExchange
{
public:
    ReplicaExchange(const double& tmin, const double& tmax,
        Molecule* original, const unsigned& Ncopies, const unsigned& sfreq,
            const int& seed);

    ~ReplicaExchange();

    void propose_moves(const unsigned& step);

    ReplicaCopy replica;

private:
    // vector of all beta's
    vector<double> beta;
    // number of copies
    unsigned Ncopies;
    // map from beta to molecule
    unordered_map<double, Molecule*> bmap;
    // uniform generator
    UniformGenerator uniform_gen;
    // switch frequency
    const unsigned sfreq;

    void propose_move(Molecule* m1, Molecule* m2);

    void linear_beta(const double& tmin, const double& tmax,
        const unsigned& N);

    bool accept_reject(const double& a);

    void update_momentum(Molecule* molecule_pt, const double& tobeta,
        const double& frombeta);

    void rescale_momentum(Molecule* molecule, const double& beta_at,
        const double& beta_new);

    void make_beta_to_mol_map();
};
#endif
