#include "observable.h"

observable::observable(const param& p_ref, const polymerSpecies& poly_ref, const ion& ion_ref)
    : p(p_ref), poly(poly_ref), ions(ion_ref), stress_obs(p_ref, poly_ref, ion_ref)
{
    build_measure_terms();
}

observable::~observable()
{
    release();
}

bool observable::initialize()
{
    release();
    build_measure_terms();
    if (!stress_obs.initialize()) return false;
    if (!allocate_isotropic_pressure_scratch()) return false;
    return true;
}

void observable::release()
{
    stress_obs.release();
    device_free(iso_block_psi);
    device_free(iso_work_a);
    device_free(iso_work_b);
    device_free(iso_work_c);
}

void observable::build_measure_terms()
{
    short_range_self_pressure_value = p.short_range_self_pressure();
}

bool observable::isotropic_pressure_enabled() const
{
    return p.stress.enabled && p.stress.isotropic_pressure;
}

bool observable::allocate_isotropic_pressure_scratch()
{
    if (!isotropic_pressure_enabled()) return true;
    const size_t size = static_cast<size_t>(p.size());
    device_alloc(iso_block_psi, size);
    device_alloc(iso_work_a, size);
    device_alloc(iso_work_b, size);
    device_alloc(iso_work_c, size);
    return true;
}

observableRecord observable::measure(cufftHandle plan,
                                     deviceState& d,
                                     const std::vector<Complex>& Qp,
                                     const Complex* omega,
                                     const Complex* psi) const
{
    observableRecord record;
    record.chemical = chemical_operators(plan, d, Qp, psi);
    record.chemical_self = self_chemical_operators();
    record.chemical_ideal = ideal_chemical_operators();
    record.hamiltonian_density = mean_field_hamiltonian_density(plan, d, record.chemical, omega, psi);
    record.stress_values = stress_obs.measure(plan, d, Qp, omega, psi);
    if (isotropic_pressure_enabled()) record.stress_values.push_back(isotropic_pressure_operator(plan, d, Qp, omega, psi));
    return record;
}

Complex observable::chemical_operator(Complex Q) const
{
    return complex_neg_log(Q);
}

Complex observable::ion_chemical(cufftHandle plan, deviceState& d, const Complex* psi, int species_index) const
{
    if (!p.electrostatics_enabled()) return complex_zero();
    const int size = p.size();
    const std::vector<IonSpecies>& list = ions.all();
    const double* kernel = d.ion_kernels + static_cast<size_t>(species_index) * size;
    return chemical_operator(ion_partition_function(plan,
                                                    d.obs_a,
                                                    d.obs_b,
                                                    psi,
                                                    kernel,
                                                    list[static_cast<size_t>(species_index)].valence,
                                                    size));
}

Complex observable::polymer_self_chemical(int species_index) const
{
    const polymerSpecies::Species& sp = poly.species(species_index);
    const double electrostatic = -p.lB * sp.charge_square_sum() / (2.0 * std::sqrt(Pi) * p.a);
    return make_complex(electrostatic, 0.0);
}

Complex observable::ion_self_chemical(int species_index) const
{
    const IonSpecies& species = ions.all()[static_cast<size_t>(species_index)];
    return make_complex(ion_self_chemical_potential(p.lB, species.valence, species.smearing_length), 0.0);
}

std::vector<Complex> observable::chemical_operators(cufftHandle plan, deviceState& d, const std::vector<Complex>& Qp, const Complex* psi) const
{
    std::vector<Complex> values;
    values.reserve(Qp.size() + ions.count());
    for (Complex Q : Qp) values.push_back(chemical_operator(Q));
    for (size_t j = 0; j < ions.count(); ++j) values.push_back(ion_chemical(plan, d, psi, static_cast<int>(j)));
    return values;
}

std::vector<Complex> observable::self_chemical_operators() const
{
    std::vector<Complex> values;
    values.reserve(static_cast<size_t>(poly.count()) + ions.count());
    for (int i = 0; i < poly.count(); ++i) values.push_back(polymer_self_chemical(i));
    for (size_t j = 0; j < ions.count(); ++j) values.push_back(ion_self_chemical(static_cast<int>(j)));
    return values;
}

std::vector<Complex> observable::ideal_chemical_operators() const
{
    std::vector<Complex> values;
    values.reserve(static_cast<size_t>(poly.count()) + ions.count());
    for (const polymerSpecies::Species& sp : poly.all()) values.push_back(make_complex(ideal_chemical_potential(sp.num_chains, p.volume()), 0.0));
    for (const IonSpecies& species : ions.all()) values.push_back(make_complex(ideal_chemical_potential(species.count, p.volume()), 0.0));
    return values;
}

double observable::short_range_self_pressure() const
{
    return short_range_self_pressure_value;
}

std::vector<std::string> observable::stress_labels() const
{
    std::vector<std::string> labels = stress_obs.labels();
    if (isotropic_pressure_enabled()) labels.push_back("pressure_iso");
    return labels;
}

Complex observable::field_square_density(const Complex* field, double coefficient, deviceState& d) const
{
    const int size = p.size();
    copy_complex(d.obs_a, field, size);
    multiply_complex(d.obs_a, field, size);
    return complex_scale(sum_complex(d.obs_a, size), coefficient / static_cast<double>(size));
}

Complex observable::electrostatic_gradient_density(cufftHandle plan, deviceState& d, const Complex* psi) const
{
    if (!p.electrostatics_enabled()) return complex_zero();
    const int size = p.size();
    copy_complex(d.obs_a, psi, size);
    fft_smooth(plan, d.obs_a, d.k2, size);
    copy_complex(d.obs_b, psi, size);
    multiply_complex(d.obs_b, d.obs_a, size);
    return complex_scale(sum_complex(d.obs_b, size), 1.0 / (8.0 * Pi * p.lB * static_cast<double>(size)));
}

Complex observable::mean_field_hamiltonian_density(cufftHandle plan,
                                                   deviceState& d,
                                                   const std::vector<Complex>& chemical,
                                                   const Complex* omega,
                                                   const Complex* psi) const
{
    Complex h = complex_zero();
    const int size = p.size();
    for (int a = 0; a < p.omega_count(); ++a) {
        h = complex_add(h, field_square_density(omega + static_cast<size_t>(a) * size, 0.5 * p.omega_coefficient(a), d));
    }
    if (p.electrostatics_enabled()) h = complex_add(h, electrostatic_gradient_density(plan, d, psi));

    size_t offset = 0;
    for (const polymerSpecies::Species& sp : poly.all()) {
        if (offset < chemical.size()) h = complex_add(h, complex_scale(chemical[offset], static_cast<double>(sp.num_chains) / p.volume()));
        ++offset;
    }
    for (size_t j = 0; j < ions.count(); ++j) {
        if (offset + j < chemical.size()) h = complex_add(h, complex_scale(chemical[offset + j], static_cast<double>(ions.all()[j].count) / p.volume()));
    }
    return h;
}

void observable::build_isotropic_trace_block(cufftHandle plan,
                                             Complex* out,
                                             const Complex* field,
                                             double smearing_length_squared,
                                             double local_coefficient,
                                             const double* k2,
                                             int size) const
{
    fft_apply_real_kernel(plan, field, out, k2, size);
    scale_complex_real_const(out, -smearing_length_squared / 3.0, size);
    add_scaled_complex(out, field, local_coefficient, size);
}

Complex observable::isotropic_ideal_self_term() const
{
    long long ion_count = 0;
    for (const IonSpecies& species : ions.all()) ion_count += species.count;
    const double ideal_pressure = (static_cast<double>(poly.total_chains()) + static_cast<double>(ion_count)) / p.volume();
    return make_complex(ideal_pressure + short_range_self_pressure(), 0.0);
}

Complex observable::isotropic_polymer_field_term(cufftHandle plan,
                                                 deviceState& d,
                                                 const polymerSpecies::Species& sp,
                                                 Complex Qp) const
{
    const int size = p.size();
    const double a2 = p.a * p.a;
    Complex total = complex_zero();

    for (int s = 0; s < sp.N; ++s) {
        const Complex* qF_s = d.qF + static_cast<size_t>(s) * size;
        const Complex* qB_s = d.qB + static_cast<size_t>(s) * size;
        const Complex* expW_s = d.exp_W + static_cast<size_t>(s) * size;
        const int alpha = sp.monomer_type(s);

        build_isotropic_trace_block(plan,
                                    iso_work_b,
                                    d.W_type_s + static_cast<size_t>(alpha) * size,
                                    a2,
                                    0.5,
                                    d.k2,
                                    size);
        if (p.electrostatics_enabled()) {
            build_isotropic_trace_block(plan, iso_work_c, d.psi_s, a2, 1.0 / 6.0, d.k2, size);
            add_scaled_i_complex(iso_work_b, iso_work_c, sp.charge(s), size);
        }

        copy_complex(iso_work_a, qF_s, size);
        multiply_complex(iso_work_a, qB_s, size);
        multiply_complex(iso_work_a, expW_s, size);
        multiply_complex(iso_work_a, iso_work_b, size);
        const Complex avg = complex_scale(sum_complex(iso_work_a, size), 1.0 / static_cast<double>(size));
        total = complex_add(total, avg);
    }

    const Complex prefactor = complex_scale(complex_inverse(Qp), static_cast<double>(sp.num_chains) / p.volume());
    return complex_mul(prefactor, total);
}

Complex observable::isotropic_chain_connectivity_term(cufftHandle plan,
                                                      deviceState& d,
                                                      const polymerSpecies::Species& sp,
                                                      Complex Qp) const
{
    const int size = p.size();
    Complex total = complex_zero();

    for (int s = 0; s < sp.N - 1; ++s) {
        const Complex* qF_s = d.qF + static_cast<size_t>(s) * size;
        const Complex* qB_s = d.qB + static_cast<size_t>(s) * size;
        const Complex* expW_s = d.exp_W + static_cast<size_t>(s) * size;

        copy_complex(iso_work_a, expW_s, size);
        multiply_complex(iso_work_a, qB_s, size);
        fft_apply_real_kernel(plan, iso_work_a, iso_work_b, d.k2, size);
        copy_complex(iso_work_c, qF_s, size);
        multiply_complex(iso_work_c, iso_work_b, size);

        const Complex avg = complex_scale(sum_complex(iso_work_c, size), 1.0 / static_cast<double>(size));
        const double b2_over_9 = sp.bond_length(s) * sp.bond_length(s) / 9.0;
        total = complex_add(total, complex_scale(avg, b2_over_9));
    }

    const Complex prefactor = complex_scale(complex_inverse(Qp), static_cast<double>(sp.num_chains) / p.volume());
    return complex_mul(prefactor, total);
}

Complex observable::isotropic_ion_field_term(cufftHandle plan,
                                             deviceState& d,
                                             const Complex* psi) const
{
    if (!p.electrostatics_enabled()) return complex_zero();
    const int size = p.size();
    Complex total = complex_zero();

    for (size_t j = 0; j < ions.count(); ++j) {
        const IonSpecies& species = ions.all()[j];
        const double* kernel = d.ion_kernels + j * static_cast<size_t>(size);
        const Complex Qj = ion_partition_function(plan, d.ion_psi, d.ion_boltzmann, psi, kernel, species.valence, size);
        const double alpha2 = species.smearing_length * species.smearing_length;
        build_isotropic_trace_block(plan, iso_work_b, d.ion_psi, alpha2, 1.0 / 6.0, d.k2, size);

        copy_complex(iso_work_c, d.ion_boltzmann, size);
        multiply_complex(iso_work_c, iso_work_b, size);
        const Complex avg = complex_scale(sum_complex(iso_work_c, size), 1.0 / static_cast<double>(size));
        Complex prefactor = complex_scale(complex_inverse(Qj), static_cast<double>(species.count) * species.valence / p.volume());
        prefactor = complex_mul(make_complex(0.0, 1.0), prefactor);
        total = complex_add(total, complex_mul(prefactor, avg));
    }
    return total;
}

Complex observable::isotropic_pressure_operator(cufftHandle plan,
                                                deviceState& d,
                                                const std::vector<Complex>& Qp,
                                                const Complex* omega,
                                                const Complex* psi) const
{
    Complex value = isotropic_ideal_self_term();

    poly.build_smoothed_sources(plan, d, p, omega, psi);
    const int species_count = poly.count();
    for (int idx = 0; idx < species_count; ++idx) {
        const Complex solved_Q = poly.solve_chain_species_prepared(plan, d, p, idx);
        const Complex Q = idx < static_cast<int>(Qp.size()) ? Qp[static_cast<size_t>(idx)] : solved_Q;
        const polymerSpecies::Species& sp = poly.species(idx);
        value = complex_add(value, isotropic_polymer_field_term(plan, d, sp, Q));
        value = complex_add(value, isotropic_chain_connectivity_term(plan, d, sp, Q));
    }

    value = complex_add(value, isotropic_ion_field_term(plan, d, psi));
    return value;
}
