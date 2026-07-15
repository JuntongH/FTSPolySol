#include "stress.h"

namespace {

bool is_diagonal(stress::Component c)
{
    return c == stress::XX || c == stress::YY || c == stress::ZZ;
}

const char* component_label(stress::Component c)
{
    switch (c) {
    case stress::XX: return "stress_xx";
    case stress::YY: return "stress_yy";
    case stress::ZZ: return "stress_zz";
    case stress::XY: return "stress_xy";
    case stress::XZ: return "stress_xz";
    case stress::YZ: return "stress_yz";
    default: return "stress_unknown";
    }
}

Complex mean_product(Complex* scratch, const Complex* a, const Complex* b, int size)
{
    copy_complex(scratch, a, size);
    multiply_complex(scratch, b, size);
    return complex_scale(sum_complex(scratch, size), 1.0 / static_cast<double>(size));
}

} // namespace

stress::stress(const param& p_ref, const polymerSpecies& poly_ref, const ion& ion_ref)
    : p(p_ref), poly(poly_ref), ions(ion_ref)
{
}

stress::~stress()
{
    release();
}

bool stress::enabled() const
{
    return p.stress.enabled && !active_components.empty();
}

size_t stress::count() const
{
    return active_components.size();
}

void stress::build_active_components()
{
    active_components.clear();
    if (!p.stress.enabled) return;
    if (p.stress.component_xx) active_components.push_back(XX);
    if (p.stress.component_yy) active_components.push_back(YY);
    if (p.stress.component_zz) active_components.push_back(ZZ);
    if (p.stress.component_xy) active_components.push_back(XY);
    if (p.stress.component_xz) active_components.push_back(XZ);
    if (p.stress.component_yz) active_components.push_back(YZ);
}

std::vector<std::string> stress::labels() const
{
    std::vector<std::string> out;
    out.reserve(count());
    for (Component c : active_components) out.push_back(component_label(c));
    return out;
}

bool stress::initialize()
{
    release();
    build_active_components();
    if (!p.stress.enabled || active_components.empty()) return true;
    if (!build_kernels()) return false;
    if (!allocate_scratch()) return false;
    return true;
}

void stress::release()
{
    device_free(kab);
    device_free(minus_kab);
    device_free(scratch_a);
    device_free(scratch_b);
    device_free(scratch_c);
    device_free(dfield);
    device_free(dpsi);
    active_components.clear();
    for (int c = 0; c < NumComponents; ++c) mode_term[c] = 0.0;
}

bool stress::allocate_scratch()
{
    const size_t size = static_cast<size_t>(p.size());
    device_alloc(scratch_a, size);
    device_alloc(scratch_b, size);
    device_alloc(scratch_c, size);
    device_alloc(dfield, size);
    device_alloc(dpsi, size);
    return true;
}

bool stress::build_kernels()
{
    const int size = p.size();
    if (size <= 0) return errHandle::invalid("stress mesh size");

    std::vector<double> h_kab(static_cast<size_t>(NumComponents) * size, 0.0);
    std::vector<double> h_minus_kab(static_cast<size_t>(NumComponents) * size, 0.0);
    double sum_kab_over_k2[NumComponents] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    long long mode_count = 0;

    for (int x = 0; x < p.m1; ++x) {
        const double kx = wave_number(x, p.m1, p.L1);
        for (int y = 0; y < p.m2; ++y) {
            const double ky = wave_number(y, p.m2, p.L2);
            for (int z = 0; z < p.m3; ++z) {
                const double kz = wave_number(z, p.m3, p.L3);
                const int i = p.index(x, y, z);
                const double values[NumComponents] = {kx * kx, ky * ky, kz * kz, kx * ky, kx * kz, ky * kz};
                const double k2 = kx * kx + ky * ky + kz * kz;
                for (int c = 0; c < NumComponents; ++c) {
                    h_kab[static_cast<size_t>(c) * size + i] = values[c];
                    h_minus_kab[static_cast<size_t>(c) * size + i] = -values[c];
                }
                if (k2 > 0.0) {
                    ++mode_count;
                    for (int c = 0; c < NumComponents; ++c) sum_kab_over_k2[c] += values[c] / k2;
                }
            }
        }
    }

    for (int c = 0; c < NumComponents; ++c) {
        const double delta = is_diagonal(static_cast<Component>(c)) ? 1.0 : 0.0;
        mode_term[c] = p.electrostatics_enabled()
                         ? (((1.0 / 3.0) * static_cast<double>(mode_count) * delta - sum_kab_over_k2[c]) / p.volume())
                         : 0.0;
    }

    device_alloc(kab, h_kab.size());
    device_alloc(minus_kab, h_minus_kab.size());
    copy_to_device(kab, h_kab);
    copy_to_device(minus_kab, h_minus_kab);
    return true;
}

const double* stress::kab_component(Component c) const
{
    return kab + static_cast<size_t>(c) * static_cast<size_t>(p.size());
}

const double* stress::minus_kab_component(Component c) const
{
    return minus_kab + static_cast<size_t>(c) * static_cast<size_t>(p.size());
}

std::vector<Complex> stress::measure(cufftHandle plan,
                                     deviceState& d,
                                     const std::vector<Complex>& Qp,
                                     const Complex* omega,
                                     const Complex* psi) const
{
    std::vector<Complex> values(active_components.size(), complex_zero());
    if (!enabled()) return values;

    for (size_t i = 0; i < active_components.size(); ++i) {
        const Component c = active_components[i];
        values[i] = ideal_and_measure_term(c);
        if (p.electrostatics_enabled()) values[i] = complex_add(values[i], electrostatic_field_term(plan, d, psi, c));
    }

    poly.build_smoothed_sources(plan, d, p, omega, psi);
    const int species_count = poly.count();
    for (int idx = 0; idx < species_count; ++idx) {
        const Complex solved_Q = poly.solve_chain_species_prepared(plan, d, p, idx);
        const Complex Q = idx < static_cast<int>(Qp.size()) ? Qp[static_cast<size_t>(idx)] : solved_Q;
        const polymerSpecies::Species& sp = poly.species(idx);
        for (size_t i = 0; i < active_components.size(); ++i) {
            const Component c = active_components[i];
            values[i] = complex_add(values[i], chain_connectivity_term(plan, d, sp, Q, c));
            values[i] = complex_add(values[i], polymer_field_coupling_term(plan, d, sp, Q, c));
        }
    }

    if (p.electrostatics_enabled()) {
        for (size_t i = 0; i < active_components.size(); ++i) {
            values[i] = complex_add(values[i], ion_field_coupling_term(plan, d, psi, active_components[i]));
        }
    }

    return values;
}

Complex stress::ideal_and_measure_term(Component c) const
{
    const double delta = is_diagonal(c) ? 1.0 : 0.0;
    long long ion_count = 0;
    for (const IonSpecies& species : ions.all()) ion_count += species.count;
    const double ideal_pressure = (static_cast<double>(poly.total_chains()) + static_cast<double>(ion_count)) / p.volume();
    const double short_range_self_pressure = p.short_range_self_pressure();
    const double value = (ideal_pressure + short_range_self_pressure) * delta + mode_term[c];
    return make_complex(value, 0.0);
}

Complex stress::electrostatic_field_term(cufftHandle plan, deviceState& d, const Complex* psi, Component c) const
{
    if (!p.electrostatics_enabled()) return complex_zero();
    const int size = p.size();
    const double delta = is_diagonal(c) ? 1.0 : 0.0;

    fft_apply_real_kernel(plan, psi, scratch_a, d.k2, size);
    const Complex grad2_mean = mean_product(scratch_b, psi, scratch_a, size);
    fft_apply_real_kernel(plan, psi, scratch_a, kab_component(c), size);
    const Complex grad_ab_mean = mean_product(scratch_b, psi, scratch_a, size);

    const Complex bracket = complex_sub(grad_ab_mean, complex_scale(grad2_mean, delta / 3.0));
    return complex_scale(bracket, 1.0 / (4.0 * Pi * p.lB));
}

Complex stress::chain_connectivity_term(cufftHandle plan,
                                        deviceState& d,
                                        const polymerSpecies::Species& sp,
                                        Complex Qp,
                                        Component c) const
{
    const int size = p.size();
    Complex total = complex_zero();

    for (int s = 0; s < sp.N - 1; ++s) {
        const Complex* qF_s = d.qF + static_cast<size_t>(s) * size;
        const Complex* qB_s = d.qB + static_cast<size_t>(s) * size;
        const Complex* expW_s = d.exp_W + static_cast<size_t>(s) * size;

        copy_complex(scratch_a, expW_s, size);
        multiply_complex(scratch_a, qB_s, size);
        fft_apply_real_kernel(plan, scratch_a, scratch_b, minus_kab_component(c), size);
        copy_complex(scratch_c, qF_s, size);
        multiply_complex(scratch_c, scratch_b, size);

        const Complex avg = complex_scale(sum_complex(scratch_c, size), 1.0 / static_cast<double>(size));
        const double b2_over_3 = sp.bond_length(s) * sp.bond_length(s) / 3.0;
        total = complex_add(total, complex_scale(avg, b2_over_3));
    }

    const Complex prefactor = complex_scale(complex_inverse(Qp), -static_cast<double>(sp.num_chains) / p.volume());
    return complex_mul(prefactor, total);
}

Complex stress::polymer_field_coupling_term(cufftHandle plan,
                                            deviceState& d,
                                            const polymerSpecies::Species& sp,
                                            Complex Qp,
                                            Component c) const
{
    const int size = p.size();
    const double a2 = p.a * p.a;
    Complex total = complex_zero();

    if (p.electrostatics_enabled()) fft_apply_real_kernel(plan, d.psi_s, dpsi, minus_kab_component(c), size);
    else clear_complex(dpsi, size);

    for (int s = 0; s < sp.N; ++s) {
        const Complex* qF_s = d.qF + static_cast<size_t>(s) * size;
        const Complex* qB_s = d.qB + static_cast<size_t>(s) * size;
        const Complex* expW_s = d.exp_W + static_cast<size_t>(s) * size;
        const int alpha = sp.monomer_type(s);
        const Complex* W_alpha = d.W_type_s + static_cast<size_t>(alpha) * size;

        fft_apply_real_kernel(plan, W_alpha, dfield, minus_kab_component(c), size);
        clear_complex(scratch_b, size);
        add_scaled_complex(scratch_b, dfield, a2, size);
        if (is_diagonal(c)) add_scaled_complex(scratch_b, W_alpha, 0.5, size);

        if (p.electrostatics_enabled()) {
            add_scaled_i_complex(scratch_b, dpsi, sp.charge(s) * a2, size);
            if (is_diagonal(c)) add_scaled_i_complex(scratch_b, d.psi_s, sp.charge(s) / 6.0, size);
        }

        copy_complex(scratch_a, qF_s, size);
        multiply_complex(scratch_a, qB_s, size);
        multiply_complex(scratch_a, expW_s, size);
        copy_complex(scratch_c, scratch_a, size);
        multiply_complex(scratch_c, scratch_b, size);

        const Complex avg = complex_scale(sum_complex(scratch_c, size), 1.0 / static_cast<double>(size));
        total = complex_add(total, avg);
    }

    const Complex prefactor = complex_scale(complex_inverse(Qp), static_cast<double>(sp.num_chains) / p.volume());
    return complex_mul(prefactor, total);
}

Complex stress::ion_field_coupling_term(cufftHandle plan, deviceState& d, const Complex* psi, Component c) const
{
    if (!p.electrostatics_enabled()) return complex_zero();
    const int size = p.size();
    Complex total = complex_zero();

    for (size_t j = 0; j < ions.count(); ++j) {
        const IonSpecies& species = ions.all()[j];
        const double* kernel = d.ion_kernels + j * static_cast<size_t>(size);
        const Complex Qj = ion_partition_function(plan, d.ion_psi, d.ion_boltzmann, psi, kernel, species.valence, size);

        fft_apply_real_kernel(plan, d.ion_psi, scratch_a, minus_kab_component(c), size);
        const double a2 = species.smearing_length * species.smearing_length;
        scale_complex_real_const(scratch_a, a2, size);
        if (is_diagonal(c)) add_scaled_complex(scratch_a, d.ion_psi, 1.0 / 6.0, size);

        copy_complex(scratch_b, d.ion_boltzmann, size);
        multiply_complex(scratch_b, scratch_a, size);
        const Complex avg = complex_scale(sum_complex(scratch_b, size), 1.0 / static_cast<double>(size));
        Complex prefactor = complex_scale(complex_inverse(Qj), static_cast<double>(species.count) * species.valence / p.volume());
        prefactor = complex_mul(make_complex(0.0, 1.0), prefactor);
        total = complex_add(total, complex_mul(prefactor, avg));
    }
    return total;
}
