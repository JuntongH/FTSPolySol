#pragma once

#include "common.h"

class deviceState {
public:
    Complex* omega = nullptr;       // active interaction fields, nOmega * grid_size
    Complex* omega_bar = nullptr;   // predictor active interaction fields, nOmega * grid_size
    Complex* omega_s = nullptr;     // Gaussian-smoothed active interaction fields, nOmega * grid_size
    Complex* drift = nullptr;       // current omega drifts, nOmega * grid_size
    Complex* drift_bar = nullptr;   // predictor omega drifts, nOmega * grid_size
    double* eta_omega = nullptr;    // real CL noises for active omega fields, nOmega * grid_size

    Complex* psi = nullptr;
    Complex* psi_bar = nullptr;
    Complex* psi_s = nullptr;
    Complex* dpsi = nullptr;
    Complex* dpsi_bar = nullptr;
    double* etaPsi = nullptr;

    Complex* W_type_s = nullptr;    // smeared monomer-type fields, M * grid_size
    Complex* rho_type = nullptr;    // smeared monomer-type densities, M * grid_size
    Complex* rho_c = nullptr;
    Complex* rho_c_poly = nullptr;
    Complex* rho_c_ion = nullptr;
    Complex* rho_c_old = nullptr;
    Complex* rho_p = nullptr;

    Complex* W = nullptr;
    Complex* exp_mW = nullptr;
    Complex* exp_W = nullptr;
    Complex* qF = nullptr;
    Complex* qB = nullptr;

    Complex* ion_psi = nullptr;
    Complex* ion_boltzmann = nullptr;
    Complex* ion_charge = nullptr;
    Complex* ion_rho = nullptr;

    Complex* obs_a = nullptr;
    Complex* obs_b = nullptr;
    Complex* obs_c = nullptr;

    double* Gamma = nullptr;
    double* PHI = nullptr;
    double* k2 = nullptr;
    double* ion_kernels = nullptr;
    int* finite_flag = nullptr;

    ~deviceState();
    void allocate(int size, int chain_slots, int bond_slots, int ion_count, int omega_count, int monomer_type_count);
    void release();
};
