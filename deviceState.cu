#include "deviceState.h"

deviceState::~deviceState()
{
    release();
}

void deviceState::allocate(int size, int chain_slots, int phi_slots, int ion_count, int omega_count, int monomer_type_count)
{
    const size_t volume = static_cast<size_t>(size);
    const size_t chain = volume * static_cast<size_t>(std::max(1, chain_slots));
    const size_t phi = volume * static_cast<size_t>(std::max(1, phi_slots));
    const size_t omega_volume = volume * static_cast<size_t>(std::max(1, omega_count));
    const size_t type_volume = volume * static_cast<size_t>(std::max(1, monomer_type_count));
    const size_t ion_volume = volume * static_cast<size_t>(ion_count);

    device_alloc(omega, omega_volume);
    device_alloc(omega_bar, omega_volume);
    device_alloc(omega_s, omega_volume);
    device_alloc(drift, omega_volume);
    device_alloc(drift_bar, omega_volume);
    device_alloc(eta_omega, omega_volume);

    device_alloc(psi, volume);
    device_alloc(psi_bar, volume);
    device_alloc(psi_s, volume);
    device_alloc(dpsi, volume);
    device_alloc(dpsi_bar, volume);
    device_alloc(etaPsi, volume);

    device_alloc(W_type_s, type_volume);
    device_alloc(rho_type, type_volume);
    device_alloc(rho_c, volume);
    device_alloc(rho_c_poly, volume);
    device_alloc(rho_c_ion, volume);
    device_alloc(rho_c_old, volume);
    device_alloc(rho_p, volume);

    device_alloc(W, chain);
    device_alloc(exp_mW, chain);
    device_alloc(exp_W, chain);
    device_alloc(qF, chain);
    device_alloc(qB, chain);

    device_alloc(obs_a, volume);
    device_alloc(obs_b, volume);
    device_alloc(obs_c, volume);

    if (ion_count > 0) {
        device_alloc(ion_psi, volume);
        device_alloc(ion_boltzmann, volume);
        device_alloc(ion_charge, volume);
        device_alloc(ion_rho, ion_volume);
        device_alloc(ion_kernels, ion_volume);
    }

    device_alloc(Gamma, volume);
    device_alloc(PHI, phi);
    device_alloc(k2, volume);
    device_alloc(finite_flag, 1);
}

void deviceState::release()
{
    device_free(omega);
    device_free(omega_bar);
    device_free(omega_s);
    device_free(drift);
    device_free(drift_bar);
    device_free(eta_omega);

    device_free(psi);
    device_free(psi_bar);
    device_free(psi_s);
    device_free(dpsi);
    device_free(dpsi_bar);
    device_free(etaPsi);

    device_free(W_type_s);
    device_free(rho_type);
    device_free(rho_c);
    device_free(rho_c_poly);
    device_free(rho_c_ion);
    device_free(rho_c_old);
    device_free(rho_p);

    device_free(W);
    device_free(exp_mW);
    device_free(exp_W);
    device_free(qF);
    device_free(qB);

    device_free(ion_psi);
    device_free(ion_boltzmann);
    device_free(ion_charge);
    device_free(ion_rho);
    device_free(obs_a);
    device_free(obs_b);
    device_free(obs_c);
    device_free(ion_kernels);
    device_free(Gamma);
    device_free(PHI);
    device_free(k2);
    device_free(finite_flag);
}
