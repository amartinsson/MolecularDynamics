#include "ReplicaExchange.hpp"

ReplicaExchange::ReplicaExchange(const double& tmin, const double& tmax,
    Molecule* original, const unsigned& Ncopies, const unsigned& sfreq, const int& seed)
        : replica(original, Ncopies), Ncopies(Ncopies), sfreq(sfreq),
            uniform_gen(0.0, 1.0, seed)
{
    // make linear map in beta
    linear_beta(tmin, tmax, Ncopies);

    // make map between molecules and temperatures
    make_beta_to_mol_map();
}

void ReplicaExchange::propose_moves(const unsigned& step)
{
    // propose to exchange all the replicas at neighbouring temperatures
    // only at specific times
    if(step >= 1 && step % sfreq == 0) {
        for(unsigned i=0; i<Ncopies-1; i++) {

            double beta1 = beta[i];
            double beta2 = beta[i+1];

            Molecule* m1 = bmap.at(beta1);
            Molecule* m2 = bmap.at(beta2);
            //
            // printf("molecule m1 beta = %f\n", beta1);
            // printf("molecule m2 beta = %f\n", beta2);
            propose_move(m1, m2);
        }
    }
}

void ReplicaExchange::propose_move(Molecule* m1, Molecule* m2) {
    // get potentials and temperatures and calculate acceptance
    double V1 = m1->potential();
    double V2 = m2->potential();

    double beta1 = m1->beta();
    double beta2 = m2->beta();

    double acc = exp(-(beta1 - beta2) * (V2 - V1));
    // double acc = exp(-(beta2 - beta1) * (V1 - V2));
    // printf("\tm1 (%f, %f): beta1 = %f, V1 =%f\n", m1->particle(0).q(0), m1->particle(0).q(1), beta1, V1);
    // printf("\tm2 (%f, %f): beta2 = %f, V2 =%f\n", m2->particle(0).q(0), m2->particle(0).q(1), beta2, V2);

    bool accepted = accept_reject(acc);

    // printf("\t switch (%f) was accepted = %i\n", acc, accepted);

    if(accepted)
    {
        // only exchange the position and temperaure
        // leave the momentum alone
        // rescale the momentum instead as described in the
        // original paper
        replica.exchange_pos_tem(m1, m2);

        // rescale the momentum of the particle
        rescale_momentum(m1, beta1, beta2);
        rescale_momentum(m2, beta2, beta1);

        // exchange temperatures
        m1->set_beta(beta2);
        m2->set_beta(beta1);
    }
}

void ReplicaExchange::rescale_momentum(Molecule* molecule,
    const double& beta_at, const double& beta_new)
{
    // printf("\t old momentum = %f, %f\n", molecule->particle(0).p(0), molecule->particle(0).p(1));
    molecule->particle(0).p *= sqrt(beta_at / beta_new);
    // molecule->particle(0).p *= sqrt(beta_new / beta_at);
    // printf("\t new momentum = %f, %f\n", molecule->particle(0).p(0), molecule->particle(0).p(1));
}

bool ReplicaExchange::accept_reject(const double& a)
{
    // return boolean
    bool accepted = false;

    // make a uniform
    double u = uniform_gen();

    // printf("\tu = %f, a= %f\n", u, a);

    if(u <= a) {
        // change boolean to accepted
        accepted = true;
    }
    else {/* do nothing */}

    // return boolean if it was accepted
    return accepted;
}

void ReplicaExchange::make_beta_to_mol_map()
{
    for(unsigned i=0; i<Ncopies; i++) {
        bmap.insert(make_pair(beta[i], replica(i)));

        replica(i)->set_beta(beta[i]);
    }
}

void ReplicaExchange::linear_beta(const double& tmin, const double& tmax,
    const unsigned& N)
{
    double betamin = 1.0 / tmax;
    double betamax = 1.0 / tmin;

    double dbeta = (betamax - betamin) / ((double)N - 1.0);

    for(unsigned i=0; i<N; i++) {
        beta.push_back(betamin + i * dbeta);
    }
}
