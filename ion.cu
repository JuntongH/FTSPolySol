#include "ion.h"

bool ion::build(const param& p, double total_polymer_charge)
{
    species.clear();
    const double polymer_charge = total_polymer_charge;
    const bool combine = p.salt.combine_identical_species;
    if (p.ions_enabled && p.counterions.enabled) {
        if (p.electrostatics_enabled() && p.counterions.valence == 0.0) return errHandle::invalid("counterion valence");
        long long n = p.counterions.number;
        if (p.counterions.auto_number) {
            if (p.electrostatics_enabled()) {
                const double needed = -polymer_charge / p.counterions.valence;
                if (needed < -1.0e-10) return errHandle::invalid("counterion sign");
                if (!integer_from_double(needed, n)) return errHandle::invalid("counterion number");
            } else {
                n = 0;
            }
        }
        if (n < 0 || p.counterions.smearing_length <= 0.0) return errHandle::invalid("counterion");
        append(p.counterions.name, p.counterions.valence, n, p.counterions.smearing_length, combine);
    }
    if (p.ions_enabled && p.salt.enabled) {
        if (p.salt.num_salt < 0 || p.salt.cation_stoich < 0 || p.salt.anion_stoich < 0) return errHandle::invalid("salt counts");
        if (p.salt.smearing_length <= 0.0) return errHandle::invalid("salt");
        if (p.electrostatics_enabled()) {
            if (p.salt.cation_valence == 0.0 || p.salt.anion_valence == 0.0) return errHandle::invalid("salt");
            const double salt_charge = p.salt.cation_stoich * p.salt.cation_valence + p.salt.anion_stoich * p.salt.anion_valence;
            if (std::fabs(salt_charge) > 1.0e-10) return errHandle::message("salt is not neutral");
        }
        append(p.salt.cation_name, p.salt.cation_valence, p.salt.num_salt * static_cast<long long>(p.salt.cation_stoich), p.salt.smearing_length, combine);
        append(p.salt.anion_name, p.salt.anion_valence, p.salt.num_salt * static_cast<long long>(p.salt.anion_stoich), p.salt.smearing_length, combine);
    }
    double total_charge = polymer_charge;
    double charge_scale = std::fabs(polymer_charge);
    for (const IonSpecies& s : species) {
        if (s.count < 0 || (p.electrostatics_enabled() && s.valence == 0.0) || s.smearing_length <= 0.0) return errHandle::invalid("ion species");
        total_charge += s.valence * static_cast<double>(s.count);
        charge_scale += std::fabs(s.valence * static_cast<double>(s.count));
    }
    if (p.electrostatics_enabled() &&
        std::fabs(total_charge) > 1.0e-8 * std::max(1.0, charge_scale)) {
        return errHandle::message("electroneutrality failed");
    }
    return true;
}

void ion::build_kernels(const param& p, const std::vector<double>& k2, std::vector<double>& ion_kernels) const
{
    build_smearing_kernels(p.size(), species, k2, ion_kernels);
}


void ion::compute_densities(cufftHandle plan, deviceState& d, const param& p, const Complex* psi) const
{
    const int size = p.size();
    for (size_t j = 0; j < species.size(); ++j) {
        const IonSpecies& s = species[j];
        const double* kernel = d.ion_kernels + j * static_cast<size_t>(size);
        Complex* rho = d.ion_rho + j * static_cast<size_t>(size);
        const double number_density = static_cast<double>(s.count) / p.volume();
        if (!p.electrostatics_enabled()) {
            fill_complex(rho, make_complex(number_density, 0.0), size);
            continue;
        }
        ion_density_from_psi(plan, rho, d.ion_psi, d.ion_boltzmann, psi, kernel,
                             s.valence, number_density, size);
        fft_apply_real_kernel_scaled(plan, rho, d.ion_charge, kernel, s.valence, size);
        add_complex(d.rho_c_ion, d.ion_charge, size);
        add_complex(d.rho_c, d.ion_charge, size);
    }
}

void ion::print(const param& p) const
{
    for (size_t j = 0; j < species.size(); ++j) {
        const IonSpecies& s = species[j];
        std::printf("ion[%zu] %s z=% 13.4e n=%lld smearing=% 13.4e density=% 13.4e\n", j, s.name.c_str(), s.valence, s.count, s.smearing_length, static_cast<double>(s.count) / p.volume());
    }
}

const std::vector<IonSpecies>& ion::all() const
{
    return species;
}

bool ion::empty() const
{
    return species.empty();
}

size_t ion::count() const
{
    return species.size();
}

bool ion::nearly_equal(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol;
}

bool ion::integer_from_double(double value, long long& out)
{
    const double rounded = std::round(value);
    const double scale = std::max(1.0, std::fabs(value));
    if (std::fabs(value - rounded) > 1.0e-8 * scale) return false;
    out = static_cast<long long>(rounded);
    return true;
}

std::string ion::clean_name(std::string name)
{
    for (char& c : name) if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') c = '_';
    return name.empty() ? "ion" : name;
}

void ion::append(const std::string& name, double valence, long long count_value, double smearing_length, bool combine)
{
    if (count_value == 0) return;
    if (combine) {
        for (IonSpecies& s : species) {
            if (nearly_equal(s.valence, valence, 1.0e-12) && nearly_equal(s.smearing_length, smearing_length, 1.0e-12)) {
                s.count += count_value;
                return;
            }
        }
    }
    IonSpecies s;
    s.name = clean_name(name);
    s.valence = valence;
    s.count = count_value;
    s.smearing_length = smearing_length;
    species.push_back(s);
}
