#include "Template.h"

#include <vector>

#include "Species.h"

using namespace std;


template <int Model>
template <int RateId>
inline void Template<Model>::monteCarloRoutine(Particles *particles, unsigned int ipart, vector<double> *Epart,
                                                      Patch *patch, Projector *Proj, const unsigned int Z,
                                                      const electricFields E, vector<double> &IonizRate_tunnel,
                                                      vector<double> &Dnom_tunnel) 
{
    double TotalIonizPot, factorJion, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    LocalFields Jion;
    const double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;
    unsigned int newZ, Zp1, k_times;
    factorJion = factorJion_0 * E.inv * E.inv;
    ran_p = patch->rand_->uniform();
    IonizRate_tunnel[Z] = Template<Model>::ionizationRate<RateId>(Z, E);

    // Total ionization potential (used to compute the ionization current)
    TotalIonizPot = 0.0;

    // k_times will give the nb of ionization events
    k_times = 0;
    Zp1 = Z + 1;

    if (Zp1 == atomic_number_)
    {
        // if ionization of the last electron: single ionization
        // -----------------------------------------------------
        if (ran_p < 1.0 - exp(-IonizRate_tunnel[Z] * dt))
        {
            TotalIonizPot += Potential[Z];
            k_times = 1;
        }
    }
    else
    {
        // else : multiple ionization can occur in one time-step
        //        partial & final ionization are decoupled (see Nuter Phys.
        //        Plasmas)
        // -------------------------------------------------------------------------

        // initialization
        Mult = 1.0;
        Dnom_tunnel[0] = 1.0;
        Pint_tunnel = exp(-IonizRate_tunnel[Z] * dt);  // cummulative prob.

        // multiple ionization loop while Pint_tunnel < ran_p and still partial
        // ionization
        while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1))
        {
            newZ = Zp1 + k_times;
            IonizRate_tunnel[newZ] = Template<Model>::ionizationRate<RateId>(newZ, E);
            D_sum = 0.0;
            P_sum = 0.0;
            Mult *= IonizRate_tunnel[Z + k_times];
            for (unsigned int i = 0; i < k_times + 1; i++)
            {
                Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                D_sum += Dnom_tunnel[i];
                P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
            }
            Dnom_tunnel[k_times + 1] -= D_sum;
            P_sum = P_sum + Dnom_tunnel[k_times + 1] * exp(-IonizRate_tunnel[newZ] * dt);
            Pint_tunnel = Pint_tunnel + P_sum * Mult;

            TotalIonizPot += Potential[Z + k_times];
            k_times++;
        }  // END while

        // final ionization (of last electron)
        if (((1.0 - Pint_tunnel) > ran_p) && (k_times == atomic_number_ - Zp1))
        {
            TotalIonizPot += Potential[atomic_number_ - 1];
            k_times++;
        }
    }  // END Multiple ionization routine

    // Compute ionization current
    if (patch->EMfields->Jx_ != NULL)
    {  // For the moment ionization current is
       // not accounted for in AM geometry
        factorJion *= TotalIonizPot;
        Jion.x = factorJion * *(E.x + ipart);
        Jion.y = factorJion * *(E.y + ipart);
        Jion.z = factorJion * *(E.z + ipart);

        Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion);
    }

    // Creation of the new electrons
    // (variable weights are used)
    // -----------------------------

    if (k_times != 0)
    {
        new_electrons.createParticle();
        int idNew = new_electrons.size() - 1;
        for (unsigned int i = 0; i < new_electrons.dimension(); i++)
        {
            new_electrons.position(i, idNew) = particles->position(i, ipart);
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
        }
        new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
        new_electrons.charge(idNew) = -1;

        if (save_ion_charge_)
        {
            ion_charge_.push_back(particles->charge(ipart));
        }

        // Increase the charge of the particle
        particles->charge(ipart) += k_times;
    }
}

