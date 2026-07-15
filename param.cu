#include "param.h"

#include <sstream>

namespace {

bool parse_polymer_section_index(const std::string& section, int& index)
{
    const std::string prefix = "polymer.";
    if (section.size() <= prefix.size()) return false;
    if (section.compare(0, prefix.size(), prefix) != 0) return false;
    const std::string text = section.substr(prefix.size());
    if (text.empty()) return false;
    char* end = nullptr;
    const long value = std::strtol(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0' || value < 0 || value > std::numeric_limits<int>::max()) return false;
    index = static_cast<int>(value);
    return true;
}

int bracket_delta(const std::string& text)
{
    int delta = 0;
    bool quoted = false;
    char quote = '\0';
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if ((c == '"' || c == '\'') && (i == 0 || text[i - 1] != '\\')) {
            if (!quoted) { quoted = true; quote = c; }
            else if (quote == c) { quoted = false; }
            continue;
        }
        if (quoted) continue;
        if (c == '[') ++delta;
        else if (c == ']') --delta;
    }
    return delta;
}

std::vector<std::string> bracket_tokens(std::string value)
{
    for (char& c : value) {
        if (c == '[' || c == ']' || c == ',' || c == ';') c = ' ';
    }
    std::vector<std::string> tokens;
    std::istringstream input(value);
    std::string token;
    while (input >> token) tokens.push_back(parse::unquote(token));
    return tokens;
}

bool parse_string_list(const std::string& value, std::vector<std::string>& out)
{
    out = bracket_tokens(value);
    return !out.empty();
}

bool parse_real_list(const std::string& value, std::vector<double>& out)
{
    const std::vector<std::string> tokens = bracket_tokens(value);
    out.clear();
    out.reserve(tokens.size());
    for (const std::string& token : tokens) {
        char* end = nullptr;
        const double x = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || !std::isfinite(x)) return false;
        out.push_back(x);
    }
    return !out.empty();
}

std::string default_monomer_name(int index)
{
    if (index >= 0 && index < 26) return std::string(1, static_cast<char>('A' + index));
    return std::string("T") + std::to_string(index);
}

} // namespace

int param::size() const
{
    return static_cast<int>(grid_size(m1, m2, m3));
}

double param::volume() const
{
    return L1 * L2 * L3;
}

bool param::electrostatics_enabled() const
{
    return lB > 0.0;
}

int param::index(int x, int y, int z) const
{
    return static_cast<int>(linear_index(x, y, z, m2, m3));
}

int param::monomer_index(const std::string& name) const
{
    return modes.monomer_index(name);
}

int param::omega_count() const
{
    return modes.active_count();
}

int param::monomer_count() const
{
    return modes.monomer_count();
}

double param::omega_coefficient(int active_index) const
{
    return rho0 * modes.active_inv_abs_eigenvalue(active_index);
}

double param::omega_mobility(int active_index) const
{
    return lambda_omega.at(static_cast<size_t>(active_index));
}

Complex param::omega_gamma(int active_index) const
{
    return modes.active_gamma(active_index);
}

double param::omega_O(int active_index, int alpha) const
{
    return modes.active_O(active_index, alpha);
}

double param::short_range_self_pressure() const
{
    if (a <= 0.0 || total_segments <= 0 || modes.K_matrix().empty()) return 0.0;
    const std::vector<double>& K = modes.K_matrix();
    const int M = monomer_type_count;
    long double weighted_diag = 0.0L;
    for (int alpha = 0; alpha < M && alpha < static_cast<int>(total_type_segments.size()); ++alpha) {
        weighted_diag += static_cast<long double>(total_type_segments[static_cast<size_t>(alpha)]) *
                         static_cast<long double>(K[static_cast<size_t>(alpha) * M + alpha]);
    }
    const long double mean_diag = weighted_diag / static_cast<long double>(total_segments);
    const long double denom = 16.0L * static_cast<long double>(Pi) * std::sqrt(static_cast<long double>(Pi)) *
                              static_cast<long double>(a) * static_cast<long double>(a) * static_cast<long double>(a);
    return static_cast<double>(mean_diag / denom);
}

bool param::legacy_ab_field_order() const
{
    return monomer_type_count == 2 && !monomer_count_set && !monomer_names_set && !chi_matrix_set;
}


param::polymer_species_input& param::ensure_polymer_species(int idx)
{
    if (idx < 0) idx = 0;
    if (polymer_species.size() <= static_cast<size_t>(idx)) polymer_species.resize(static_cast<size_t>(idx + 1));
    polymer_species[static_cast<size_t>(idx)].section_seen = true;
    return polymer_species[static_cast<size_t>(idx)];
}

bool param::read(const char* filename)
{
    FILE* fp = std::fopen(filename, "r");
    if (!fp) return errHandle::file(filename);
    std::string section;
    char line[4096];
    int line_number = 0;
    while (std::fgets(line, sizeof(line), fp)) {
        ++line_number;
        std::string text = parse::trim(parse::strip_comment(line));
        if (text.empty()) continue;
        if (text.front() == '[' && text.back() == ']') {
            section = parse::trim(text.substr(1, text.size() - 2));
            continue;
        }
        const size_t eq = text.find('=');
        if (eq == std::string::npos) { std::fclose(fp); return errHandle::line(filename, line_number); }

        const std::string key = parse::trim(text.substr(0, eq));
        std::string value = parse::trim(text.substr(eq + 1));
        int balance = bracket_delta(value);
        while (balance > 0 && std::fgets(line, sizeof(line), fp)) {
            ++line_number;
            const std::string extra = parse::trim(parse::strip_comment(line));
            if (!extra.empty()) {
                value += " ";
                value += extra;
                balance += bracket_delta(extra);
            }
        }
        if (balance != 0) {
            std::fclose(fp);
            return errHandle::message(std::string(filename) + ": unbalanced brackets while reading [" + section + "] " + key);
        }
        if (!apply(section, key, value)) { std::fclose(fp); return false; }
    }
    std::fclose(fp);
    return derive();
}

bool param::apply(const std::string& section, const std::string& key, const std::string& value)
{
    bool ok = true;
    int polymer_index = -1;
    if (section == "grid") {
        if (key == "m1") ok = parse::integer(value, m1);
        else if (key == "m2") ok = parse::integer(value, m2);
        else if (key == "m3") ok = parse::integer(value, m3);
        else if (key == "a") ok = parse::real(value, a);
        else if (key == "dL") ok = parse::real(value, dL);
        else return unknown(section, key);
    } else if (section == "monomers") {
        if (key == "M" || key == "m" || key == "count") { ok = parse::integer(value, monomer_type_count); monomer_count_set = ok; }
        else if (key == "names" || key == "types") { ok = parse_string_list(value, monomer_names); monomer_names_set = ok; }
        else return unknown(section, key);
    } else if (section == "polymer") {
        if (key == "species_count") { ok = parse::integer(value, species_count); species_count_set = ok; }
        else if (key == "sequence_file") ok = parse::string_value(value, sequence_file);
        else return unknown(section, key);
    } else if (parse_polymer_section_index(section, polymer_index)) {
        polymer_species_input& sp = ensure_polymer_species(polymer_index);
        if (key == "name") { ok = parse::string_value(value, sp.name); sp.name_set = ok; }
        else if (key == "num_chains") { ok = parse::integer64(value, sp.num_chains); sp.num_chains_set = ok; }
        else if (key == "N") { ok = parse::integer(value, sp.N); sp.N_set = ok; }
        else return unknown(section, key);
    } else if (section == "interactions") {
        if (key == "chi" || key == "chi_matrix") { ok = parse_real_list(value, chi_matrix); chi_matrix_set = ok; }
        else if (key == "xi" || key == "v") ok = parse::real(value, xi);
        else if (key == "lB") ok = parse::real(value, lB);
        else if (key == "eigen_tol_rel" || key == "interaction_eigen_tol_rel") ok = parse::real(value, eigen_tol_rel);
        else return unknown(section, key);
    } else if (section == "dynamics") {
        if (key == "dt") ok = parse::real(value, dt);
        else if (key == "lambda_omega" || key == "lambda_w" || key == "lambda") ok = parse_real_list(value, lambda_omega_input);
        else if (key == "lambdaPsi") ok = parse::real(value, lambdaPsi);
        else if (key == "Nt") ok = parse::integer(value, Nt);
        else if (key == "Nss") ok = parse::integer(value, Nss);
        else if (key == "Ns") ok = parse::integer(value, Ns);
        else if (key == "Ntmp" || key == "ntmp") ok = parse::integer(value, Ntmp);
        else if (key == "Nso" || key == "nso") ok = parse::integer(value, Nso);
        else if (key == "N_stat" || key == "n_stat" || key == "Nstat" || key == "nstat") ok = parse::integer(value, N_stat);
        else if (key == "max_nan_restarts") ok = parse::integer(value, max_nan_restarts);
        else if (key == "noise") ok = parse::boolean(value, noise);
        else return unknown(section, key);
    } else if (section == "ions") {
        if (key == "enabled") ok = parse::boolean(value, ions_enabled);
        else if (key == "default_smearing_length") ok = parse::real(value, default_smearing_length);
        else return unknown(section, key);
    } else if (section == "ions.counterions") {
        if (key == "enabled") ok = parse::boolean(value, counterions.enabled);
        else if (key == "number") ok = parse::counterion_number(value, counterions.auto_number, counterions.number);
        else if (key == "valence") ok = parse::real(value, counterions.valence);
        else if (key == "name") ok = parse::string_value(value, counterions.name);
        else if (key == "smearing_length") ok = parse::real(value, counterions.smearing_length);
        else return unknown(section, key);
    } else if (section == "ions.salt") {
        if (key == "enabled") ok = parse::boolean(value, salt.enabled);
        else if (key == "num_salt") ok = parse::integer64(value, salt.num_salt);
        else if (key == "cation_name") ok = parse::string_value(value, salt.cation_name);
        else if (key == "cation_valence") ok = parse::real(value, salt.cation_valence);
        else if (key == "cation_stoich") ok = parse::integer(value, salt.cation_stoich);
        else if (key == "anion_name") ok = parse::string_value(value, salt.anion_name);
        else if (key == "anion_valence") ok = parse::real(value, salt.anion_valence);
        else if (key == "anion_stoich") ok = parse::integer(value, salt.anion_stoich);
        else if (key == "smearing_length") ok = parse::real(value, salt.smearing_length);
        else if (key == "combine_identical_species") ok = parse::boolean(value, salt.combine_identical_species);
        else return unknown(section, key);
    } else if (section == "structure_factor") {
        // Compatibility only: structure-factor output has been removed from the
        // generalized interaction-mode code path.  These legacy keys are
        // accepted and ignored so older input files remain readable.
        bool bool_sink = false;
        int int_sink = 0;
        double real_sink = 0.0;
        std::string string_sink;
        if (key == "enabled") ok = parse::boolean(value, bool_sink);
        else if (key == "file" || key == "output_file") ok = parse::string_value(value, string_sink);
        else if (key == "desmear" || key == "particle_center") ok = parse::boolean(value, bool_sink);
        else if (key == "skip_zero") ok = parse::boolean(value, bool_sink);
        else if (key == "interval") ok = parse::integer(value, int_sink);
        else if (key == "kmax") ok = parse::real(value, real_sink);
        else if (key == "shell_width" || key == "bin_width" || key == "dk") ok = parse::real(value, real_sink);
        else if (key == "shell_tolerance" || key == "k_tolerance") ok = parse::real(value, real_sink);
        else return unknown(section, key);
    } else if (section == "observables") {
        if (key == "stress") ok = parse::boolean(value, stress.enabled);
        else if (key == "component_xx") ok = parse::boolean(value, stress.component_xx);
        else if (key == "component_yy") ok = parse::boolean(value, stress.component_yy);
        else if (key == "component_zz") ok = parse::boolean(value, stress.component_zz);
        else if (key == "component_xy") ok = parse::boolean(value, stress.component_xy);
        else if (key == "component_xz") ok = parse::boolean(value, stress.component_xz);
        else if (key == "component_yz") ok = parse::boolean(value, stress.component_yz);
        else if (key == "isotropic_pressure" || key == "pressure" || key == "pressure_iso" || key == "p_iso" || key == "pressure_trace" || key == "p_trace") ok = parse::boolean(value, stress.isotropic_pressure);
        else return unknown(section, key);
    } else if (section == "output") {
        bool bool_sink = false;
        double real_sink = 0.0;
        std::string string_sink;
        if (key == "Nso" || key == "nso") ok = parse::integer(value, Nso);
        else if (key == "N_stat" || key == "n_stat" || key == "Nstat" || key == "nstat") ok = parse::integer(value, N_stat);
        else if (key == "omega_file" || key == "w_file") ok = parse::string_value(value, omega_file);
        else if (key == "psi_file") ok = parse::string_value(value, psi_file);
        else if (key == "tmp_omega_file" || key == "tmp_w_file") ok = parse::string_value(value, tmp_omega_file);
        else if (key == "tmp_psi_file") ok = parse::string_value(value, tmp_psi_file);
        else if (key == "rho_type_file") ok = parse::string_value(value, rho_type_file);
        else if (key == "rhop_file") ok = parse::string_value(value, rhop_file);
        else if (key == "rhoc_file") ok = parse::string_value(value, rhoc_file);
        else if (key == "rhoc_poly_file") ok = parse::string_value(value, rhoc_poly_file);
        else if (key == "rhoc_ion_file") ok = parse::string_value(value, rhoc_ion_file);
        else if (key == "ion_density_prefix") ok = parse::string_value(value, ion_density_prefix);
        else if (key == "step_omega_file" || key == "step_w_file") ok = parse::string_value(value, step_omega_file);
        else if (key == "step_psi_file") ok = parse::string_value(value, step_psi_file);
        else if (key == "step_rho_type_file") ok = parse::string_value(value, step_rho_type_file);
        else if (key == "step_rhop_file") ok = parse::string_value(value, step_rhop_file);
        else if (key == "step_rhoc_file") ok = parse::string_value(value, step_rhoc_file);
        else if (key == "step_rhoc_poly_file") ok = parse::string_value(value, step_rhoc_poly_file);
        else if (key == "step_rhoc_ion_file") ok = parse::string_value(value, step_rhoc_ion_file);
        else if (key == "step_ion_density_prefix") ok = parse::string_value(value, step_ion_density_prefix);
        else if (key == "observable_file") ok = parse::string_value(value, observable_file);
        else if (key == "drift_max_file") ok = parse::string_value(value, drift_max_file);
        else if (key == "order_param" || key == "order_param_file") ok = parse::string_value(value, string_sink);
        else if (key == "order_param_kc" || key == "order_kc") ok = parse::real(value, real_sink);
        else if (key == "order_param_enabled" || key == "order_parameter_enabled") ok = parse::boolean(value, bool_sink);
        else return unknown(section, key);
    } else {
        return unknown(section, key);
    }
    if (!ok) return errHandle::parse("[" + section + "] " + key);
    return true;
}

bool param::unknown(const std::string& section, const std::string& key) const
{
    if (section.empty()) return errHandle::unknown("parameter " + key);
    return errHandle::unknown("parameter [" + section + "] " + key);
}

bool param::derive()
{
    if (species_count > 0 && polymer_species.size() < static_cast<size_t>(species_count)) {
        polymer_species.resize(static_cast<size_t>(species_count));
    }

    L1 = dL * static_cast<double>(m1);
    L2 = dL * static_cast<double>(m2);
    L3 = dL * static_cast<double>(m3);
    dV = dL * dL * dL;
    if (Ns <= 0) Ns = Nt;
    if (default_smearing_length == UnsetSmearing) default_smearing_length = a;
    if (counterions.smearing_length == UnsetSmearing) counterions.smearing_length = default_smearing_length;
    if (salt.smearing_length == UnsetSmearing) salt.smearing_length = default_smearing_length;
    recompute_polymer_derived();
    if (!finalize_interactions()) return false;
    return finalize_mobilities();
}

bool param::finalize_interactions()
{
    if (!monomer_count_set && monomer_names_set) {
        monomer_type_count = static_cast<int>(monomer_names.size());
    } else if (!monomer_count_set && !monomer_names_set && chi_matrix_set) {
        const int root = static_cast<int>(std::sqrt(static_cast<double>(chi_matrix.size())) + 0.5);
        if (root > 0 && static_cast<size_t>(root) * static_cast<size_t>(root) == chi_matrix.size()) monomer_type_count = root;
    } else if (!monomer_count_set && !monomer_names_set) {
        // Backward compatibility: old A/B input files did not have a
        // [monomers] section.  In that case default to the legacy A/B labels.
        monomer_type_count = 2;
    }

    if (monomer_type_count <= 0) return errHandle::invalid("[monomers] M");
    if (monomer_names_set) {
        if (monomer_names.size() != static_cast<size_t>(monomer_type_count)) return errHandle::invalid("[monomers] names length must equal M");
    } else {
        monomer_names.clear();
        for (int i = 0; i < monomer_type_count; ++i) monomer_names.push_back(default_monomer_name(i));
    }


    const size_t mm = static_cast<size_t>(monomer_type_count) * static_cast<size_t>(monomer_type_count);
    if (!chi_matrix_set) chi_matrix.assign(mm, 0.0);
    if (chi_matrix.size() != mm) return errHandle::invalid("[interactions] chi matrix size must be M*M");

    return modes.build(monomer_type_count, monomer_names, xi, chi_matrix, eigen_tol_rel);
}

bool param::finalize_mobilities()
{
    const int n = omega_count();
    lambda_omega.clear();
    if (n <= 0) return true;
    if (lambda_omega_input.empty()) {
        lambda_omega.assign(static_cast<size_t>(n), 0.0);
        return true;
    }
    if (!lambda_omega_input.empty() && lambda_omega_input.size() == 1) {
        lambda_omega.assign(static_cast<size_t>(n), lambda_omega_input[0]);
    } else if (!lambda_omega_input.empty() && lambda_omega_input.size() == static_cast<size_t>(n)) {
        lambda_omega = lambda_omega_input;
    } else if (!lambda_omega_input.empty()) {
        return errHandle::invalid("lambda_omega length must be 1 or the number of active interaction modes");
    }
    for (double x : lambda_omega) {
        if (!std::isfinite(x) || x < 0.0) return errHandle::invalid("lambda_omega");
    }
    return true;
}

void param::recompute_polymer_derived()
{
    const int n_species = std::max(0, species_count);
    if (n_species > 0 && polymer_species.size() < static_cast<size_t>(n_species)) polymer_species.resize(static_cast<size_t>(n_species));

    Nmax = 0;
    total_polymer_chains = 0;
    total_segments = 0;
    total_type_segments.assign(static_cast<size_t>(std::max(0, monomer_type_count)), 0);

    for (int idx = 0; idx < n_species; ++idx) {
        polymer_species_input& sp = polymer_species[static_cast<size_t>(idx)];
        if (sp.type_counts.size() != static_cast<size_t>(std::max(0, monomer_type_count))) {
            sp.type_counts.assign(static_cast<size_t>(std::max(0, monomer_type_count)), 0);
        }
        if (sp.N > Nmax) Nmax = sp.N;
        total_polymer_chains += sp.num_chains;
        total_segments += sp.num_chains * static_cast<long long>(sp.N);
        for (int alpha = 0; alpha < monomer_type_count; ++alpha) {
            total_type_segments[static_cast<size_t>(alpha)] += sp.num_chains * static_cast<long long>(sp.type_counts[static_cast<size_t>(alpha)]);
        }
    }

    if (volume() > 0.0 && total_segments > 0) rho0 = static_cast<double>(total_segments) / volume();
}

bool param::validate_modes() const
{
    bool ok = true;
    if (!modes.valid()) { errHandle::invalid("interaction modes"); ok = false; }
    if (monomer_type_count != modes.monomer_count()) { errHandle::invalid("monomer type count"); ok = false; }
    if (monomer_names != modes.monomer_names()) { errHandle::invalid("monomer names"); ok = false; }
    if (lambda_omega.size() != static_cast<size_t>(std::max(0, omega_count()))) { errHandle::invalid("lambda_omega count"); ok = false; }
    return ok;
}

bool param::validate() const
{
    bool ok = true;
    auto check_step_template = [&](const std::string& pattern, const char* name) {
        if (pattern.empty()) return;
        if (!valid_step_filename_template(pattern)) { errHandle::invalid(std::string(name) + " (use exactly one %d in the file name)"); ok = false; }
    };
    const long long mesh_size = static_cast<long long>(grid_size(m1, m2, m3));
    if (!validate_modes()) ok = false;
    if (!species_count_set) { errHandle::invalid("[polymer] species_count is required"); ok = false; }
    if (species_count <= 0) { errHandle::invalid("polymer species_count"); ok = false; }
    if (sequence_file.empty()) { errHandle::invalid("[polymer] sequence_file is required"); ok = false; }
    if (polymer_species.size() < static_cast<size_t>(std::max(0, species_count))) { errHandle::invalid("polymer species list"); ok = false; }
    if (polymer_species.size() > static_cast<size_t>(std::max(0, species_count))) { errHandle::invalid("polymer species section index exceeds species_count"); ok = false; }
    for (int idx = 0; idx < species_count && idx < static_cast<int>(polymer_species.size()); ++idx) {
        const polymer_species_input& sp = polymer_species[static_cast<size_t>(idx)];
        const std::string prefix = "[polymer." + std::to_string(idx) + "]";
        if (!sp.section_seen) { errHandle::invalid(prefix + " section is required"); ok = false; }
        if (!sp.name_set) { errHandle::invalid(prefix + " name is required"); ok = false; }
        if (!sp.N_set) { errHandle::invalid(prefix + " N is required"); ok = false; }
        if (!sp.num_chains_set) { errHandle::invalid(prefix + " num_chains is required"); ok = false; }
        if (sp.N <= 0) { errHandle::invalid(prefix + " N"); ok = false; }
        int count_sum = 0;
        for (int c : sp.type_counts) count_sum += c;
        if (count_sum != sp.N) { errHandle::invalid(prefix + " derived monomer type counts"); ok = false; }
        if (sp.num_chains <= 0) { errHandle::invalid(prefix + " num_chains"); ok = false; }
    }
    if (Nmax <= 0) { errHandle::invalid("Nmax"); ok = false; }
    if (total_polymer_chains <= 0) { errHandle::invalid("total polymer chains"); ok = false; }
    if (total_segments <= 0) { errHandle::invalid("total polymer segments"); ok = false; }
    if (m1 <= 0 || m2 <= 0 || m3 <= 0) { errHandle::invalid("mesh"); ok = false; }
    if (mesh_size <= 0 || mesh_size > 2147483647LL) { errHandle::invalid("mesh size"); ok = false; }
    if (a <= 0.0 || dL <= 0.0 || lB < 0.0) { errHandle::invalid("a, dL, or lB"); ok = false; }
    if (dt < 0.0 || lambdaPsi < 0.0) { errHandle::invalid("dt or lambdaPsi"); ok = false; }
    if (Nss < 0 || Ns < 0 || Ntmp < 0 || Nso < -1 || N_stat < -1) { errHandle::invalid("sampling/checkpoint/step-output/drift-stat intervals"); ok = false; }
    if (N_stat >= 0 && drift_max_file.empty()) { errHandle::invalid("drift_max_file"); ok = false; }
    if (stress.enabled && !stress.any_observable()) { errHandle::invalid("stress observable requires at least one component_* = true or isotropic_pressure = true"); ok = false; }
    if (!stress.enabled && stress.any_observable()) { errHandle::invalid("stress components or isotropic_pressure require [observables] stress = true"); ok = false; }
    check_step_template(step_omega_file, "step_omega_file");
    check_step_template(step_psi_file, "step_psi_file");
    check_step_template(step_rho_type_file, "step_rho_type_file");
    check_step_template(step_rhop_file, "step_rhop_file");
    check_step_template(step_rhoc_file, "step_rhoc_file");
    check_step_template(step_rhoc_poly_file, "step_rhoc_poly_file");
    check_step_template(step_rhoc_ion_file, "step_rhoc_ion_file");
    check_step_template(step_ion_density_prefix, "step_ion_density_prefix");
    if (max_nan_restarts < 0) { errHandle::invalid("max_nan_restarts"); ok = false; }
    if (rho0 <= 0.0) { errHandle::invalid("rho0"); ok = false; }
    if (xi < 0.0) { errHandle::invalid("xi"); ok = false; }
    if (default_smearing_length <= 0.0) { errHandle::invalid("default_smearing_length"); ok = false; }
    if (counterions.smearing_length <= 0.0) { errHandle::invalid("counterion smearing_length"); ok = false; }
    if (salt.smearing_length <= 0.0) { errHandle::invalid("salt smearing_length"); ok = false; }
    return ok;
}

void param::print() const
{
    modes.print_mode_table(stdout);
    std::printf("species_count = %d\n", species_count);
    for (int idx = 0; idx < species_count && idx < static_cast<int>(polymer_species.size()); ++idx) {
        const polymer_species_input& sp = polymer_species[static_cast<size_t>(idx)];
        std::printf("polymer.%d name = %s\n", idx, sp.name.c_str());
        std::printf("polymer.%d N = %d\n", idx, sp.N);
        for (int alpha = 0; alpha < monomer_type_count && alpha < static_cast<int>(sp.type_counts.size()); ++alpha) {
            std::printf("polymer.%d count_%s = %d\n", idx, monomer_names[static_cast<size_t>(alpha)].c_str(), sp.type_counts[static_cast<size_t>(alpha)]);
        }
        std::printf("polymer.%d num_chains = %lld\n", idx, sp.num_chains);
    }
    std::printf("Nmax = %d\n", Nmax);
    std::printf("total_polymer_chains = %lld\n", total_polymer_chains);
    std::printf("total_segments = %lld\n", total_segments);
    for (int alpha = 0; alpha < monomer_type_count && alpha < static_cast<int>(total_type_segments.size()); ++alpha) {
        std::printf("total_%s_segments = %lld\n", monomer_names[static_cast<size_t>(alpha)].c_str(), total_type_segments[static_cast<size_t>(alpha)]);
    }
    std::printf("m1 m2 m3 = %d %d %d\n", m1, m2, m3);
    std::printf("a = % 13.4e\n", a);
    std::printf("dL = % 13.4e\n", dL);
    std::printf("L1 L2 L3 = % 13.4e % 13.4e % 13.4e\n", L1, L2, L3);
    std::printf("dV = % 13.4e\n", dV);
    std::printf("V = % 13.4e\n", volume());
    std::printf("reduced density total_segments/V = % 13.4e\n", rho0);
    std::printf("xi = % 13.4e\n", xi);
    std::printf("short_range_self_pressure = % 13.4e\n", short_range_self_pressure());
    std::printf("lB = % 13.4e%s\n", lB, electrostatics_enabled() ? "" : " (electrostatics disabled)");
    for (int a = 0; a < omega_count(); ++a) {
        std::printf("omega_mode[%d] coefficient rho0/|d| = % 13.4e lambda = % 13.4e\n", a, omega_coefficient(a), omega_mobility(a));
    }
    std::printf("dt = % 13.4e\n", dt);
    std::printf("lambdaPsi = % 13.4e\n", lambdaPsi);
    std::printf("Nt = %d\n", Nt);
    std::printf("Nss = %d\n", Nss);
    std::printf("Ns = %d\n", Ns);
    std::printf("Ntmp = %d\n", Ntmp);
    std::printf("Nso = %d\n", Nso);
    if (N_stat >= 0) {
        std::printf("N_stat = %d\n", N_stat);
        std::printf("drift_max_file = %s\n", drift_max_file.c_str());
    }
    if (stress.enabled) {
        std::printf("stress = on\n");
        std::printf("stress_component_xx = %s\n", stress.component_xx ? "on" : "off");
        std::printf("stress_component_yy = %s\n", stress.component_yy ? "on" : "off");
        std::printf("stress_component_zz = %s\n", stress.component_zz ? "on" : "off");
        std::printf("stress_component_xy = %s\n", stress.component_xy ? "on" : "off");
        std::printf("stress_component_xz = %s\n", stress.component_xz ? "on" : "off");
        std::printf("stress_component_yz = %s\n", stress.component_yz ? "on" : "off");
        std::printf("stress_isotropic_pressure = %s\n", stress.isotropic_pressure ? "on" : "off");
    }
    if (!step_omega_file.empty()) std::printf("step_omega_file = %s\n", step_omega_file.c_str());
    if (!step_psi_file.empty()) std::printf("step_psi_file = %s\n", step_psi_file.c_str());
    if (!step_rho_type_file.empty()) std::printf("step_rho_type_file = %s\n", step_rho_type_file.c_str());
    if (!step_rhop_file.empty()) std::printf("step_rhop_file = %s\n", step_rhop_file.c_str());
    if (!step_rhoc_file.empty()) std::printf("step_rhoc_file = %s\n", step_rhoc_file.c_str());
    if (!step_rhoc_poly_file.empty()) std::printf("step_rhoc_poly_file = %s\n", step_rhoc_poly_file.c_str());
    if (!step_rhoc_ion_file.empty()) std::printf("step_rhoc_ion_file = %s\n", step_rhoc_ion_file.c_str());
    if (!step_ion_density_prefix.empty()) std::printf("step_ion_density_prefix = %s\n", step_ion_density_prefix.c_str());
    std::printf("max_nan_restarts = %d\n", max_nan_restarts);
    std::printf("noise = %s\n", noise ? "on" : "off");
    std::printf("ions = %s\n", ions_enabled ? "on" : "off");
}
