#pragma once

#include "polymer.h"
#include "ion.h"
#include "stress.h"

struct observableRecord {
    std::vector<Complex> chemical;
    std::vector<Complex> chemical_self;
    std::vector<Complex> chemical_ideal;
    std::vector<Complex> stress_values;
    Complex hamiltonian_density;
};

class observable {
public:
    observable(const param& p, const polymerSpecies& poly, const ion& ions);
    ~observable();

    bool initialize();
    void release();

    observableRecord measure(cufftHandle plan,
                             deviceState& d,
                             const std::vector<Complex>& Qp,
                             const Complex* omega,
                             const Complex* psi) const;
    double short_range_self_pressure() const;
    std::vector<std::string> stress_labels() const;

private:
    const param& p;
    const polymerSpecies& poly;
    const ion& ions;
    stress stress_obs;
    double short_range_self_pressure_value = 0.0;

    mutable Complex* iso_block_psi = nullptr;
    mutable Complex* iso_work_a = nullptr;
    mutable Complex* iso_work_b = nullptr;
    mutable Complex* iso_work_c = nullptr;

    void build_measure_terms();
    bool isotropic_pressure_enabled() const;
    bool allocate_isotropic_pressure_scratch();
    Complex chemical_operator(Complex Q) const;
    Complex ion_chemical(cufftHandle plan, deviceState& d, const Complex* psi, int species_index) const;
    Complex polymer_self_chemical(int species_index) const;
    Complex ion_self_chemical(int species_index) const;
    std::vector<Complex> chemical_operators(cufftHandle plan, deviceState& d, const std::vector<Complex>& Qp, const Complex* psi) const;
    std::vector<Complex> self_chemical_operators() const;
    std::vector<Complex> ideal_chemical_operators() const;
    Complex mean_field_hamiltonian_density(cufftHandle plan, deviceState& d, const std::vector<Complex>& chemical, const Complex* omega, const Complex* psi) const;
    Complex field_square_density(const Complex* field, double coefficient, deviceState& d) const;
    Complex electrostatic_gradient_density(cufftHandle plan, deviceState& d, const Complex* psi) const;
    void build_isotropic_trace_block(cufftHandle plan,
                                     Complex* out,
                                     const Complex* field,
                                     double smearing_length_squared,
                                     double local_coefficient,
                                     const double* k2,
                                     int size) const;
    Complex isotropic_pressure_operator(cufftHandle plan,
                                        deviceState& d,
                                        const std::vector<Complex>& Qp,
                                        const Complex* omega,
                                        const Complex* psi) const;
    Complex isotropic_ideal_self_term() const;
    Complex isotropic_polymer_field_term(cufftHandle plan,
                                         deviceState& d,
                                         const polymerSpecies::Species& sp,
                                         Complex Qp) const;
    Complex isotropic_chain_connectivity_term(cufftHandle plan,
                                              deviceState& d,
                                              const polymerSpecies::Species& sp,
                                              Complex Qp) const;
    Complex isotropic_ion_field_term(cufftHandle plan,
                                     deviceState& d,
                                     const Complex* psi) const;
};
