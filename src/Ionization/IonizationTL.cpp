#include "IonizationTunnel.h"
#include "IonizationTables.h"

template <>
IonizationTunnel<2>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tong-Lin Tunnel Ionizaton class");

    ionization_tl_parameter_ =
        species->ionization_tl_parameter_;  // species->ionization_tl_parameter_ is double
                                            // Varies from 6 to 9. This is the alpha parameter in Tong-Lin exponential, see Eq. (6) in [M F Ciappina and S V Popruzhenko 2020 Laser Phys. Lett. 17 025301 2020].
    lambda_tunnel.resize(atomic_number_);

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (int Z = 0; Z < (int)atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);
        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
        lambda_tunnel[Z] = ionization_tl_parameter_ * pow(cst, 2) / gamma_tunnel[Z];
    }

    DEBUG("Finished Creating the Tong-Lin Tunnel Ionizaton class");
}

template <>
double IonizationTunnel<2>::ionizationRate(const int Z, const double E, const int oldZ)
{
    const double delta = gamma_tunnel[Z] / E;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E * lambda_tunnel[Z]);
}
