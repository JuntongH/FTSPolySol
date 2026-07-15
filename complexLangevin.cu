#include "complexLangevin.h"
#include "fieldIO.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr int StepColumnWidth = 10;
constexpr int ValueColumnWidth = 22;
constexpr int ValuePrecision = 10;
constexpr int ScreenValueWidth = 13;
constexpr int ScreenPrecision = 4;

void append_ion_labels(std::vector<std::string>& labels, const ion& ions, const char* suffix)
{
    const std::vector<IonSpecies>& list = ions.all();
    for (size_t j = 0; j < list.size(); ++j) {
        labels.push_back("mu_ion" + std::to_string(j) + "_" + list[j].name + "_" + suffix);
    }
}

std::string polymer_label(const param& p, int species_index, const char* suffix)
{
    if (p.species_count == 1) return std::string("mu_poly_") + suffix;
    std::string name = species_index < static_cast<int>(p.polymer_species.size())
                           ? p.polymer_species[static_cast<size_t>(species_index)].name
                           : ("P" + std::to_string(species_index));
    for (char& c : name) if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') c = '_';
    return "mu_poly" + std::to_string(species_index) + "_" + name + "_" + suffix;
}

void append_polymer_labels(std::vector<std::string>& labels, const param& p, const char* suffix)
{
    for (int i = 0; i < p.species_count; ++i) labels.push_back(polymer_label(p, i, suffix));
}

std::vector<std::string> observable_labels(const param& p, const ion& ions, const observable& obs)
{
    const std::vector<std::string> stress_labels = obs.stress_labels();
    std::vector<std::string> labels;
    labels.reserve(3 * static_cast<size_t>(p.species_count) + 3 * ions.count() + stress_labels.size() + 2);
    append_polymer_labels(labels, p, "ex");
    append_polymer_labels(labels, p, "id");
    append_polymer_labels(labels, p, "self");
    append_ion_labels(labels, ions, "ex");
    append_ion_labels(labels, ions, "id");
    append_ion_labels(labels, ions, "self");
    labels.push_back("Hmf_density_real");
    labels.push_back("Hmf_density_imag");
    labels.insert(labels.end(), stress_labels.begin(), stress_labels.end());
    return labels;
}

std::vector<int> observable_widths(const std::vector<std::string>& labels)
{
    std::vector<int> widths;
    widths.reserve(labels.size());
    for (const std::string& label : labels) widths.push_back(std::max(ValueColumnWidth, static_cast<int>(label.size())));
    return widths;
}

size_t observable_column_count(const param& p, const ion& ions, const observable& obs)
{
    return observable_labels(p, ions, obs).size();
}

std::vector<double> observable_values(const observableRecord& record, int polymer_count)
{
    const size_t polymer_values = static_cast<size_t>(std::max(0, polymer_count));
    std::vector<double> values;
    values.reserve(record.chemical.size() + record.chemical_self.size() + record.chemical_ideal.size() + record.stress_values.size() + 2);

    for (size_t i = 0; i < polymer_values && i < record.chemical.size(); ++i) values.push_back(record.chemical[i].x);
    for (size_t i = 0; i < polymer_values && i < record.chemical_ideal.size(); ++i) values.push_back(record.chemical_ideal[i].x);
    for (size_t i = 0; i < polymer_values && i < record.chemical_self.size(); ++i) values.push_back(record.chemical_self[i].x);

    for (size_t j = polymer_values; j < record.chemical.size(); ++j) values.push_back(record.chemical[j].x);
    for (size_t j = polymer_values; j < record.chemical_ideal.size(); ++j) values.push_back(record.chemical_ideal[j].x);
    for (size_t j = polymer_values; j < record.chemical_self.size(); ++j) values.push_back(record.chemical_self[j].x);

    values.push_back(record.hamiltonian_density.x);
    values.push_back(record.hamiltonian_density.y);

    for (const Complex& value : record.stress_values) values.push_back(value.x);
    return values;
}

void write_observable_labels(FILE* fp, const param& p, const ion& ions, const observable& obs, const char* first_label)
{
    const std::vector<std::string> labels = observable_labels(p, ions, obs);
    const std::vector<int> widths = observable_widths(labels);
    std::fprintf(fp, "#%*s", StepColumnWidth - 1, first_label);
    for (size_t i = 0; i < labels.size(); ++i) std::fprintf(fp, " %*s", widths[i], labels[i].c_str());
}

void write_value(FILE* fp, double value, int width)
{
    std::fprintf(fp, " % *.*e", width, ValuePrecision, value);
}

void write_observable_row(FILE* fp, const param& p, const ion& ions, const observable& obs, int first_value, const std::vector<double>& values)
{
    const std::vector<std::string> labels = observable_labels(p, ions, obs);
    const std::vector<int> widths = observable_widths(labels);
    std::fprintf(fp, "%*d", StepColumnWidth, first_value);
    for (size_t i = 0; i < values.size(); ++i) write_value(fp, values[i], i < widths.size() ? widths[i] : ValueColumnWidth);
}

void write_observable_row(FILE* fp, const param& p, const ion& ions, const observable& obs, const char* first_value, const std::vector<double>& values)
{
    const std::vector<std::string> labels = observable_labels(p, ions, obs);
    const std::vector<int> widths = observable_widths(labels);
    std::fprintf(fp, "%*s", StepColumnWidth, first_value);
    for (size_t i = 0; i < values.size(); ++i) write_value(fp, values[i], i < widths.size() ? widths[i] : ValueColumnWidth);
}

void write_screen_value(double value)
{
    std::printf("% *.*e", ScreenValueWidth, ScreenPrecision, value);
}

bool finite_observable_record(const observableRecord& record, int polymer_count, int step)
{
    for (size_t i = 0; i < record.chemical.size(); ++i) {
        const std::string name = "chemical[" + std::to_string(i) + "]";
        if (!finite_complex_value(name.c_str(), record.chemical[i], step)) return false;
    }
    for (size_t i = 0; i < record.chemical_self.size(); ++i) {
        const std::string name = "chemical_self[" + std::to_string(i) + "]";
        if (!finite_complex_value(name.c_str(), record.chemical_self[i], step)) return false;
    }
    for (size_t i = 0; i < record.chemical_ideal.size(); ++i) {
        const std::string name = "chemical_ideal[" + std::to_string(i) + "]";
        if (!finite_complex_value(name.c_str(), record.chemical_ideal[i], step)) return false;
    }
    for (size_t i = 0; i < record.stress_values.size(); ++i) {
        const std::string name = "stress[" + std::to_string(i) + "]";
        if (!finite_complex_value(name.c_str(), record.stress_values[i], step)) return false;
    }
    const std::vector<double> values = observable_values(record, polymer_count);
    for (size_t i = 0; i < values.size(); ++i) {
        const std::string name = "observable_value[" + std::to_string(i) + "]";
        if (!finite_real_value(name.c_str(), values[i], step)) return false;
    }
    if (!finite_complex_value("hamiltonian_density", record.hamiltonian_density, step)) return false;
    return true;
}

} // namespace

complexLangevin::complexLangevin(param& param_ref, polymerSpecies& polymer_ref, ion& ion_ref)
    : p(param_ref), poly(polymer_ref), ions(ion_ref), obs(param_ref, polymer_ref, ion_ref)
{
}

complexLangevin::~complexLangevin()
{
    release();
}

bool complexLangevin::initialize(const std::vector<Complex>& omega,
                                 const std::vector<Complex>& psi,
                                 const std::vector<double>& Gamma,
                                 const std::vector<double>& PHI,
                                 const std::vector<double>& k2,
                                 const std::vector<double>& ion_kernels)
{
    const size_t omega_expected = static_cast<size_t>(std::max(0, p.omega_count())) * static_cast<size_t>(p.size());
    if (omega.size() != omega_expected) return errHandle::message("omega input size does not match active interaction modes");
    if (psi.size() != static_cast<size_t>(p.size())) return errHandle::message("psi input size does not match mesh");

    d.allocate(p.size(), std::max(1, p.Nmax), static_cast<int>(poly.bond_kernel_count()), static_cast<int>(ions.count()), p.omega_count(), p.monomer_count());
    fft_plan_3d(&plan, p.m1, p.m2, p.m3);
    ready = true;

    copy_to_device(d.omega, omega);
    copy_to_device(d.psi, psi);
    remove_psi_zero_mode(d.psi);
    copy_to_device(d.Gamma, Gamma);
    copy_to_device(d.PHI, PHI);
    copy_to_device(d.k2, k2);
    if (!ions.empty()) copy_to_device(d.ion_kernels, ion_kernels);

    if (!obs.initialize()) return false;
    return true;
}

bool complexLangevin::run(bool screenprint)
{
    if (!open_observables()) return false;
    if (!open_drift_statistics()) return false;

    int target_step = p.Nt;
    if (checkpointing_enabled()) {
        bool nonfinite = false;
        if (!write_checkpoint(0, nonfinite)) return false;
    }

    int step = 0;
    while (step < target_step) {
        const int next_step = step + 1;
        const StepStatus status = run_step(screenprint, next_step);
        if (status == StepStatus::failed) return false;
        if (status == StepStatus::nonfinite) {
            if (!restart_from_checkpoint(next_step, target_step)) return false;
            step = checkpoint_step;
            continue;
        }

        step = next_step;
        if (checkpointing_enabled() && step % p.Ntmp == 0) {
            bool nonfinite = false;
            if (!write_checkpoint(step, nonfinite)) {
                if (!nonfinite) return false;
                if (!restart_from_checkpoint(step, target_step)) return false;
                step = checkpoint_step;
                continue;
            }
        }
    }

    return write_observable_average();
}

bool complexLangevin::final_densities()
{
    const std::vector<Complex> Qp = compute_densities(d.omega, d.psi);
    return finite_chain_partitions(Qp, "Q(final)", -1);
}

void complexLangevin::copy_results(std::vector<Complex>& omega,
                                   std::vector<Complex>& psi,
                                   std::vector<Complex>& rho_type,
                                   std::vector<Complex>& rho_p,
                                   std::vector<Complex>& rho_c,
                                   std::vector<Complex>& rho_c_poly,
                                   std::vector<Complex>& rho_c_ion,
                                   std::vector<Complex>& ion_rho)
{
    copy_to_host(omega, d.omega);
    copy_to_host(psi, d.psi);
    copy_to_host(rho_type, d.rho_type);
    copy_to_host(rho_p, d.rho_p);
    copy_to_host(rho_c, d.rho_c);
    copy_to_host(rho_c_poly, d.rho_c_poly);
    copy_to_host(rho_c_ion, d.rho_c_ion);
    if (!ions.empty()) copy_to_host(ion_rho, d.ion_rho);
}

std::vector<Complex> complexLangevin::compute_densities(const Complex* omega, const Complex* psi)
{
    std::vector<Complex> Qp = poly.solve_and_compute_densities(plan, d, p, omega, psi);
    ions.compute_densities(plan, d, p, psi);
    return Qp;
}

bool complexLangevin::finite_chain_partitions(const std::vector<Complex>& Qp, const char* label, int step) const
{
    if (Qp.size() != static_cast<size_t>(poly.count())) return errHandle::message(std::string(label) + " polymer species count mismatch");
    for (size_t i = 0; i < Qp.size(); ++i) {
        const std::string name = std::string(label) + "[" + std::to_string(i) + "]";
        if (!finite_complex_value(name.c_str(), Qp[i], step)) return false;
    }
    return true;
}

void complexLangevin::compute_drifts(const Complex* omega_term, Complex* drift)
{
    const int size = p.size();
    for (int a = 0; a < p.omega_count(); ++a) {
        Complex* u = drift + static_cast<size_t>(a) * size;
        const Complex* w = omega_term + static_cast<size_t>(a) * size;
        clear_complex(u, size);
        add_scaled_complex(u, w, p.omega_coefficient(a), size);
        for (int alpha = 0; alpha < p.monomer_count(); ++alpha) {
            const Complex coeff = complex_scale(p.omega_gamma(a), p.omega_O(a, alpha));
            add_scaled_complex_const(u, d.rho_type + static_cast<size_t>(alpha) * size, coeff, size);
        }
    }
}

void complexLangevin::generate_noises()
{
    const int size = p.size();
    for (int a = 0; a < p.omega_count(); ++a) {
        const double sigma = langevin_noise_stddev(p.noise, p.omega_mobility(a), p.dt, p.dV);
        generate_noise(d.eta_omega + static_cast<size_t>(a) * size, sigma, size);
    }
    const double sigmaPsi = p.electrostatics_enabled() ? langevin_noise_stddev(p.noise, p.lambdaPsi, p.dt, p.dV) : 0.0;
    generate_noise(d.etaPsi, sigmaPsi, size);
    if (p.electrostatics_enabled()) remove_real_mean(d.etaPsi, size);
}

void complexLangevin::compute_psi_drift(const Complex* psi_field, const Complex* rho_source, Complex* drift)
{
    const int size = p.size();
    if (!p.electrostatics_enabled()) {
        clear_complex(drift, size);
        return;
    }
    copy_complex(drift, psi_field, size);
    fft_smooth(plan, drift, d.k2, size);
    scale_complex_real_const(drift, 1.0 / (4.0 * Pi * p.lB), size);
    add_i_times_complex(drift, rho_source, size);
    remove_psi_zero_mode(drift);
}

void complexLangevin::remove_psi_zero_mode(Complex* psi_field)
{
    if (p.electrostatics_enabled()) remove_complex_mean(psi_field, p.size());
    else clear_complex(psi_field, p.size());
}


double complexLangevin::corrector_drift_max(Complex* scratch, const Complex* first, const Complex* second, double mobility)
{
    const int size = p.size();
    copy_complex(scratch, first, size);
    add_complex(scratch, second, size);
    return max_complex_magnitude(scratch, 0.5 * mobility, size);
}

void complexLangevin::predictor_step()
{
    const int size = p.size();
    compute_psi_drift(d.psi, d.rho_c_old, d.dpsi);

    for (int a = 0; a < p.omega_count(); ++a) {
        Complex* pred = d.omega_bar + static_cast<size_t>(a) * size;
        const Complex* old = d.omega + static_cast<size_t>(a) * size;
        const Complex* u = d.drift + static_cast<size_t>(a) * size;
        copy_complex(pred, old, size);
        add_scaled_complex(pred, u, -p.dt * p.omega_mobility(a), size);
        add_scaled_real(pred, d.eta_omega + static_cast<size_t>(a) * size, 1.0, size);
    }

    copy_complex(d.psi_bar, d.psi, size);
    if (p.electrostatics_enabled()) add_scaled_complex(d.psi_bar, d.dpsi, -p.dt * p.lambdaPsi, size);
    else clear_complex(d.psi_bar, size);
    if (p.electrostatics_enabled()) add_scaled_real(d.psi_bar, d.etaPsi, 1.0, size);
    remove_psi_zero_mode(d.psi_bar);
}

void complexLangevin::corrector_step()
{
    const int size = p.size();
    for (int a = 0; a < p.omega_count(); ++a) {
        Complex* w = d.omega + static_cast<size_t>(a) * size;
        const Complex* u0 = d.drift + static_cast<size_t>(a) * size;
        const Complex* u1 = d.drift_bar + static_cast<size_t>(a) * size;
        const double lambda = p.omega_mobility(a);
        add_scaled_complex(w, u0, -0.5 * p.dt * lambda, size);
        add_scaled_complex(w, u1, -0.5 * p.dt * lambda, size);
        add_scaled_real(w, d.eta_omega + static_cast<size_t>(a) * size, 1.0, size);
    }

    if (p.electrostatics_enabled()) {
        add_scaled_complex(d.psi, d.dpsi, -0.5 * p.dt * p.lambdaPsi, size);
        add_scaled_complex(d.psi, d.dpsi_bar, -0.5 * p.dt * p.lambdaPsi, size);
        add_scaled_real(d.psi, d.etaPsi, 1.0, size);
    }
    remove_psi_zero_mode(d.psi);
}

bool complexLangevin::open_observables()
{
    if (p.observable_file.empty()) return true;
    observable_output = std::fopen(p.observable_file.c_str(), "w");
    if (!observable_output) return errHandle::file(p.observable_file);
    observable_sum.assign(observable_column_count(p, ions, obs), 0.0);
    observable_count = 0;
    write_observable_header();
    return true;
}

bool complexLangevin::open_drift_statistics()
{
    if (!drift_statistics_enabled()) return true;
    drift_output = std::fopen(p.drift_max_file.c_str(), "w");
    if (!drift_output) return errHandle::file(p.drift_max_file);
    write_drift_header();
    return true;
}

void complexLangevin::write_observable_header()
{
    if (!observable_output) return;
    std::fprintf(observable_output, "# observables Nss=%d Ns=%d\n", p.Nss, p.Ns);
    write_observable_labels(observable_output, p, ions, obs, "step");
    std::fprintf(observable_output, "\n");
    std::fflush(observable_output);
}

void complexLangevin::write_drift_header()
{
    if (!drift_output) return;
    std::fprintf(drift_output, "#%*s", StepColumnWidth - 1, "step");
    for (int a = 0; a < p.omega_count(); ++a) {
        const std::string label = "u_omega" + std::to_string(a) + "_max";
        std::fprintf(drift_output, " %*s", ValueColumnWidth, label.c_str());
    }
    std::fprintf(drift_output, " %*s\n", ValueColumnWidth, "u_psi_max");
    std::fflush(drift_output);
}

void complexLangevin::write_observable_record(int step, const observableRecord& record)
{
    if (!observable_output) return;
    const std::vector<double> values = observable_values(record, p.species_count);
    if (values.size() != observable_sum.size()) {
        errHandle::message("observable output column size mismatch");
        return;
    }
    write_observable_row(observable_output, p, ions, obs, step, values);
    std::fprintf(observable_output, "\n");
    std::fflush(observable_output);
}

bool complexLangevin::accumulate_observable_record(int step, const observableRecord& record)
{
    if (step <= p.Nss || observable_count >= p.Ns || observable_sum.empty()) return true;

    const std::vector<std::string> labels = observable_labels(p, ions, obs);
    const std::vector<double> values = observable_values(record, p.species_count);
    if (labels.size() != values.size() || observable_sum.size() != values.size()) return errHandle::message("observable accumulator size mismatch");

    for (size_t column = 0; column < values.size(); ++column) {
        const std::string name = labels[column] + "_sum";
        if (!finite_real_value(name.c_str(), values[column], step)) return false;
        observable_sum[column] += values[column];
        if (!finite_real_value(name.c_str(), observable_sum[column], step)) return false;
    }

    ++observable_count;
    return true;
}

bool complexLangevin::write_observable_average()
{
    if (!observable_output) return true;
    std::fprintf(observable_output, "\n# ensemble_average Nss=%d Ns=%d count=%lld\n", p.Nss, p.Ns, observable_count);
    if (observable_count == 0) {
        std::fprintf(observable_output, "# no_observables\n");
        std::fflush(observable_output);
        return true;
    }

    std::vector<double> averages;
    averages.reserve(observable_sum.size());
    for (size_t i = 0; i < observable_sum.size(); ++i) {
        const std::string name = "observable_sum[" + std::to_string(i) + "]";
        if (!finite_real_value(name.c_str(), observable_sum[i], -1)) return false;
        averages.push_back(observable_sum[i] / static_cast<double>(observable_count));
    }

    write_observable_labels(observable_output, p, ions, obs, "stat");
    std::fprintf(observable_output, "\n");
    write_observable_row(observable_output, p, ions, obs, "average", averages);
    std::fprintf(observable_output, "\n");
    std::fflush(observable_output);
    return true;
}

bool complexLangevin::drift_statistics_enabled() const
{
    return p.N_stat >= 0 && !p.drift_max_file.empty();
}

bool complexLangevin::write_drift_max(int step)
{
    if (!drift_output || step <= p.N_stat) return true;

    std::vector<double> values;
    values.reserve(static_cast<size_t>(p.omega_count()) + 1);
    for (int a = 0; a < p.omega_count(); ++a) {
        const double u = corrector_drift_max(d.obs_a,
                                             d.drift + static_cast<size_t>(a) * p.size(),
                                             d.drift_bar + static_cast<size_t>(a) * p.size(),
                                             p.omega_mobility(a));
        const std::string name = "u_omega" + std::to_string(a) + "_max";
        if (!finite_real_value(name.c_str(), u, step)) return false;
        values.push_back(u);
    }
    const double u_psi = corrector_drift_max(d.obs_a, d.dpsi, d.dpsi_bar, p.lambdaPsi);
    if (!finite_real_value("u_psi_max", u_psi, step)) return false;
    values.push_back(u_psi);

    std::fprintf(drift_output, "%*d", StepColumnWidth, step);
    for (double value : values) write_value(drift_output, value, ValueColumnWidth);
    std::fprintf(drift_output, "\n");
    std::fflush(drift_output);
    return true;
}

bool complexLangevin::reset_drift_statistics_after_restart(int failed_step)
{
    if (drift_output) {
        std::fclose(drift_output);
        drift_output = nullptr;
    }
    if (!drift_statistics_enabled()) return true;

    drift_output = std::fopen(p.drift_max_file.c_str(), "w");
    if (!drift_output) return errHandle::file(p.drift_max_file);
    write_drift_header();
    std::fprintf(drift_output,
                 "# restarted after non-finite value at step %d from checkpoint step %d; previous drift statistics discarded\n",
                 failed_step, checkpoint_step);
    std::fflush(drift_output);
    return true;
}

bool complexLangevin::checkpointing_enabled() const
{
    return p.Ntmp > 0;
}

bool complexLangevin::write_step_outputs(int step, bool& nonfinite)
{
    nonfinite = false;
    if (!fieldIO::step_outputs_enabled(p)) return true;
    const int first_output_step = p.Nso > 1 ? p.Nso : 1;
    if (step < first_output_step) return true;
    const int output_index = step - first_output_step + 1;

    const bool want_omega = !p.step_omega_file.empty();
    const bool want_psi = !p.step_psi_file.empty();
    const bool want_rho_type = !p.step_rho_type_file.empty();
    const bool want_rhop = !p.step_rhop_file.empty();
    const bool want_rhoc = !p.step_rhoc_file.empty();
    const bool want_rhoc_poly = !p.step_rhoc_poly_file.empty();
    const bool want_rhoc_ion = !p.step_rhoc_ion_file.empty();
    const bool want_ion_rho = !p.step_ion_density_prefix.empty();
    const bool want_density = want_rho_type || want_rhop || want_rhoc || want_rhoc_poly || want_rhoc_ion || want_ion_rho;

    if (want_density) {
        const std::vector<Complex> Qp = compute_densities(d.omega, d.psi);
        if (!finite_chain_partitions(Qp, "Q(step output)", step)) {
            nonfinite = true;
            return false;
        }
    }

    const int size = p.size();
    std::vector<Complex> h_omega;
    std::vector<Complex> h_psi;
    std::vector<Complex> h_rho_type;
    std::vector<Complex> h_rho_p;
    std::vector<Complex> h_rho_c;
    std::vector<Complex> h_rho_c_poly;
    std::vector<Complex> h_rho_c_ion;
    std::vector<Complex> h_ion_rho;

    if (want_omega) { h_omega.resize(static_cast<size_t>(p.omega_count()) * size); copy_to_host(h_omega, d.omega); }
    if (want_psi) { h_psi.resize(size); copy_to_host(h_psi, d.psi); }
    if (want_rho_type) { h_rho_type.resize(static_cast<size_t>(p.monomer_count()) * size); copy_to_host(h_rho_type, d.rho_type); }
    if (want_rhop) { h_rho_p.resize(size); copy_to_host(h_rho_p, d.rho_p); }
    if (want_rhoc) { h_rho_c.resize(size); copy_to_host(h_rho_c, d.rho_c); }
    if (want_rhoc_poly) { h_rho_c_poly.resize(size); copy_to_host(h_rho_c_poly, d.rho_c_poly); }
    if (want_rhoc_ion) { h_rho_c_ion.resize(size); copy_to_host(h_rho_c_ion, d.rho_c_ion); }
    if (want_ion_rho && !ions.empty()) { h_ion_rho.resize(static_cast<size_t>(size) * ions.count()); copy_to_host(h_ion_rho, d.ion_rho); }

    if (!finite_complex_vectors({{"step-output omega", &h_omega},
                                {"step-output psi", &h_psi},
                                {"step-output rho_type", &h_rho_type},
                                {"step-output rho_p", &h_rho_p},
                                {"step-output rho_c", &h_rho_c},
                                {"step-output rho_c_poly", &h_rho_c_poly},
                                {"step-output rho_c_ion", &h_rho_c_ion},
                                {"step-output ion_rho", &h_ion_rho}}, step)) {
        nonfinite = true;
        return false;
    }

    return fieldIO::write_step_outputs(p, ions, output_index, h_omega, h_psi, h_rho_type, h_rho_p, h_rho_c, h_rho_c_poly, h_rho_c_ion, h_ion_rho);
}

bool complexLangevin::write_checkpoint(int step, bool& nonfinite)
{
    nonfinite = false;
    const int size = p.size();
    std::vector<Complex> h_omega(static_cast<size_t>(p.omega_count()) * size);
    std::vector<Complex> h_psi(size);
    copy_to_host(h_omega, d.omega);
    copy_to_host(h_psi, d.psi);

    if (!finite_complex_vectors({{"checkpoint omega", &h_omega}, {"checkpoint psi", &h_psi}}, step)) {
        nonfinite = true;
        return false;
    }

    const std::string new_omega = p.tmp_omega_file + ".new";
    const std::string new_psi = p.tmp_psi_file + ".new";
    if (!fieldIO::write_omega_fields(p, new_omega, h_omega)) return false;
    if (!fieldIO::write_complex_field(p, new_psi, h_psi)) {
        std::remove(new_omega.c_str());
        return false;
    }
    if (!replace_file(new_omega, p.tmp_omega_file)) {
        std::remove(new_psi.c_str());
        return false;
    }
    if (!replace_file(new_psi, p.tmp_psi_file)) return false;
    checkpoint_valid = true;
    checkpoint_step = step;
    return true;
}

bool complexLangevin::read_checkpoint(int failed_step)
{
    std::vector<Complex> h_omega;
    std::vector<Complex> h_psi;
    if (!fieldIO::read_omega_fields_exact(p, p.tmp_omega_file, h_omega)) return false;
    if (!fieldIO::read_complex_field_exact(p, p.tmp_psi_file, h_psi)) return false;
    copy_to_device(d.omega, h_omega);
    copy_to_device(d.psi, h_psi);
    remove_psi_zero_mode(d.psi);
    return device_fields_finite("checkpoint restart", failed_step);
}

bool complexLangevin::restart_from_checkpoint(int failed_step, int& target_step)
{
    if (!checkpointing_enabled() || !checkpoint_valid || checkpoint_step < 0) return errHandle::message("non-finite value encountered before a valid checkpoint was written");
    if (restart_count >= p.max_nan_restarts) return errHandle::message("maximum non-finite restarts exceeded");
    ++restart_count;
    if (!read_checkpoint(failed_step)) return false;
    if (!reset_observables_after_restart(failed_step)) return false;
    if (!reset_drift_statistics_after_restart(failed_step)) return false;
    target_step = p.Nt;
    return true;
}

bool complexLangevin::reset_observables_after_restart(int failed_step)
{
    if (observable_output) {
        std::fclose(observable_output);
        observable_output = nullptr;
    }
    observable_sum.clear();
    observable_count = 0;
    if (p.observable_file.empty()) return true;
    observable_output = std::fopen(p.observable_file.c_str(), "w");
    if (!observable_output) return errHandle::file(p.observable_file);
    observable_sum.assign(observable_column_count(p, ions, obs), 0.0);
    write_observable_header();
    std::fprintf(observable_output,
                 "# restarted after non-finite value at step %d from checkpoint step %d; previous observable records discarded\n",
                 failed_step, checkpoint_step);
    std::fflush(observable_output);
    return true;
}

bool complexLangevin::device_fields_finite(const char* label, int step)
{
    const bool omega_ok = finite_device_field_array(d.omega, p.omega_count(), p.size(), d.finite_flag);
    const bool psi_ok = finite_device_field_array(d.psi, 1, p.size(), d.finite_flag);
    if (!omega_ok || !psi_ok) {
        if (step >= 0) return errHandle::message(std::string(label) + " has non-finite value at step " + std::to_string(step));
        return errHandle::message(std::string(label) + " has non-finite value");
    }
    return true;
}

complexLangevin::StepStatus complexLangevin::run_step(bool screenprint, int step)
{
    std::vector<Complex> Qp = compute_densities(d.omega, d.psi);
    if (!finite_chain_partitions(Qp, "Q", step)) return StepStatus::nonfinite;
    if (!device_fields_finite("density fields", step)) return StepStatus::nonfinite;
    copy_complex(d.rho_c_old, d.rho_c, p.size());

    generate_noises();
    compute_drifts(d.omega, d.drift);
    predictor_step();

    Qp = compute_densities(d.omega_bar, d.psi_bar);
    if (!finite_chain_partitions(Qp, "Qbar", step)) return StepStatus::nonfinite;
    compute_psi_drift(d.psi_bar, d.rho_c, d.dpsi_bar);

    observableRecord record = obs.measure(plan, d, Qp, d.omega_bar, d.psi_bar);
    if (!finite_observable_record(record, p.species_count, step)) return StepStatus::nonfinite;
    if (!accumulate_observable_record(step, record)) return StepStatus::failed;
    if (step % std::max(1, p.Nss) == 0) write_observable_record(step, record);
    if (screenprint) print_observables(record);

    compute_drifts(d.omega_bar, d.drift_bar);
    if (!write_drift_max(step)) return StepStatus::nonfinite;
    corrector_step();

    bool nonfinite = false;
    if (!write_step_outputs(step, nonfinite)) return nonfinite ? StepStatus::nonfinite : StepStatus::failed;
    if (!device_fields_finite("updated fields", step)) return StepStatus::nonfinite;
    return StepStatus::ok;
}

void complexLangevin::print_observables(const observableRecord& record)
{
    const std::vector<double> values = observable_values(record, p.species_count);
    for (double value : values) write_screen_value(value);
    std::printf("\n");
}

void complexLangevin::release()
{
    if (observable_output) {
        std::fclose(observable_output);
        observable_output = nullptr;
    }
    if (drift_output) {
        std::fclose(drift_output);
        drift_output = nullptr;
    }
    obs.release();
    d.release();
    if (ready) {
        fft_destroy(plan);
        ready = false;
        plan = 0;
    }
}
