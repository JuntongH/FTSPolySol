#include "polymer.h"

#include <sstream>

namespace {

bool is_separator(char c)
{
    return std::isspace(static_cast<unsigned char>(c)) || c == ',' || c == ';' || c == '[' || c == ']';
}

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

std::vector<std::string> split_tokens(const std::string& text)
{
    std::string clean = text;
    for (char& c : clean) if (is_separator(c)) c = ' ';
    std::istringstream input(clean);
    std::vector<std::string> tokens;
    std::string token;
    while (input >> token) tokens.push_back(parse::unquote(token));
    return tokens;
}

bool all_names_single_character(const param& p)
{
    for (const std::string& name : p.monomer_names) if (name.size() != 1) return false;
    return true;
}

bool parse_type_line(const std::string& line, const param& p, std::vector<int>& types)
{
    types.clear();
    std::vector<std::string> tokens = split_tokens(line);
    if (tokens.size() == 1 && tokens[0].size() > 1 && all_names_single_character(p)) {
        const std::string compact = tokens[0];
        tokens.clear();
        for (char c : compact) tokens.push_back(std::string(1, c));
    }
    for (const std::string& token : tokens) {
        const int idx = p.monomer_index(token);
        if (idx < 0) return errHandle::invalid("unknown monomer label in sequence: " + token);
        types.push_back(idx);
    }
    if (types.empty()) return errHandle::message("sequence has no monomer line");
    return true;
}

void skip_charge_separators(const std::string& text, size_t& i)
{
    while (i < text.size() && is_separator(text[i])) ++i;
}

bool parse_charge_value(const std::string& text, size_t& i, double& charge)
{
    const char c = text[i];
    if (c == '+' || c == '-') {
        const double sign = (c == '+') ? 1.0 : -1.0;
        ++i;
        if (i >= text.size() || is_separator(text[i]) || text[i] == '+' || text[i] == '-') {
            charge = sign;
            return true;
        }
        const char* begin = text.c_str() + i;
        char* end = nullptr;
        const double magnitude = std::strtod(begin, &end);
        if (end == begin) return false;
        charge = sign * magnitude;
        i = static_cast<size_t>(end - text.c_str());
        return true;
    }
    const char* begin = text.c_str() + i;
    char* end = nullptr;
    charge = std::strtod(begin, &end);
    if (end == begin) return false;
    i = static_cast<size_t>(end - text.c_str());
    return true;
}

bool parse_real_text(const std::string& text, std::vector<double>& values)
{
    values.clear();
    size_t i = 0;
    while (true) {
        skip_charge_separators(text, i);
        if (i >= text.size()) break;
        double value = 0.0;
        if (!parse_charge_value(text, i, value)) return errHandle::parse("sequence real value");
        values.push_back(value);
    }
    return true;
}

bool parse_key_value_line(const std::string& text, std::string& key, std::string& value)
{
    const size_t eq = text.find('=');
    if (eq == std::string::npos) return false;
    key = parse::trim(text.substr(0, eq));
    value = parse::trim(text.substr(eq + 1));
    for (char& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return !key.empty();
}

struct SequenceBlock {
    bool present = false;
    std::string monomers;
    std::string charges;
    std::string bonds;
};

bool build_species_from_sequence(const param& p,
                                 const param::polymer_species_input& input,
                                 const SequenceBlock& block,
                                 polymerSpecies::Species& out)
{
    std::vector<int> types;
    std::vector<double> charges;
    std::vector<double> bonds;

    if (!parse_type_line(block.monomers, p, types)) return false;
    if (!parse_real_text(block.charges, charges)) return false;
    if (!block.bonds.empty() && !parse_real_text(block.bonds, bonds)) return false;

    if (input.name.empty()) return errHandle::message("polymer species name is required before reading sequence");
    if (input.num_chains <= 0) return errHandle::message("num_chains must be positive for " + input.name);
    const int N = input.N;
    if (N <= 0) return errHandle::message("N must be positive for " + input.name);
    if (static_cast<int>(types.size()) != N) return errHandle::message("sequence length mismatch for " + input.name);
    if (static_cast<int>(charges.size()) != N) return errHandle::message("charge sequence length mismatch for " + input.name);
    if (!bonds.empty() && static_cast<int>(bonds.size()) != std::max(0, N - 1)) return errHandle::message("bond sequence length must be N-1 for " + input.name);

    out.name = input.name;
    out.num_chains = input.num_chains;
    out.N = N;
    out.type_counts.assign(static_cast<size_t>(p.monomer_count()), 0);
    out.segments.clear();
    out.segments.reserve(static_cast<size_t>(N));
    for (int s = 0; s < N; ++s) {
        double bond = 1.0;
        if (s < N - 1 && !bonds.empty()) bond = bonds[static_cast<size_t>(s)];
        if (bond <= 0.0) return errHandle::invalid("sequence bond length");
        const int type = types[static_cast<size_t>(s)];
        if (type < 0 || type >= p.monomer_count()) return errHandle::invalid("sequence monomer type");
        ++out.type_counts[static_cast<size_t>(type)];
        out.segments.push_back({type, charges[static_cast<size_t>(s)], bond});
    }
    return true;
}

} // namespace

int polymerSpecies::Species::charged_segment_count() const
{
    int count = 0;
    for (const Segment& segment : segments) if (segment.charge != 0.0) ++count;
    return count;
}

double polymerSpecies::Species::charge_square_sum() const
{
    double total = 0.0;
    for (const Segment& segment : segments) total += segment.charge * segment.charge;
    return total;
}

double polymerSpecies::Species::net_chain_charge() const
{
    double total = 0.0;
    for (const Segment& segment : segments) total += segment.charge;
    return total;
}

bool polymerSpecies::setup(param& p)
{
    species_list.clear();
    total_bonds = 0;
    if (p.sequence_file.empty()) return errHandle::message("[polymer] sequence_file is required");
    if (!read_block_sequence_file(p.sequence_file.c_str(), p)) return false;
    update_counts(p);
    assign_bond_offsets();
    return true;
}

int polymerSpecies::count() const { return static_cast<int>(species_list.size()); }
const polymerSpecies::Species& polymerSpecies::species(int idx) const { return species_list[static_cast<size_t>(idx)]; }
const std::vector<polymerSpecies::Species>& polymerSpecies::all() const { return species_list; }

long long polymerSpecies::total_chains() const
{
    long long total = 0;
    for (const Species& sp : species_list) total += sp.num_chains;
    return total;
}


size_t polymerSpecies::bond_kernel_count() const { return total_bonds; }

double polymerSpecies::total_polymer_charge() const
{
    double total = 0.0;
    for (const Species& sp : species_list) total += static_cast<double>(sp.num_chains) * sp.net_chain_charge();
    return total;
}



bool polymerSpecies::read_block_sequence_file(const char* filename, param& p)
{
    FILE* fp = std::fopen(filename, "r");
    if (!fp) return errHandle::file(filename);

    std::vector<SequenceBlock> blocks(static_cast<size_t>(std::max(1, p.species_count)));
    int current = -1;
    char line[4096];
    int line_number = 0;
    while (std::fgets(line, sizeof(line), fp)) {
        ++line_number;
        const std::string clean = parse::trim(parse::strip_comment(line));
        if (clean.empty()) continue;
        if (clean.front() == '[' && clean.back() == ']') {
            const std::string section = parse::trim(clean.substr(1, clean.size() - 2));
            int idx = -1;
            if (!parse_polymer_section_index(section, idx)) { std::fclose(fp); return errHandle::unknown("sequence section [" + section + "]"); }
            if (idx < 0 || idx >= p.species_count) { std::fclose(fp); return errHandle::invalid("sequence polymer species index"); }
            current = idx;
            blocks[static_cast<size_t>(current)].present = true;
            continue;
        }
        if (current < 0) { std::fclose(fp); return errHandle::line(filename, line_number); }
        std::string key;
        std::string value;
        if (!parse_key_value_line(clean, key, value)) { std::fclose(fp); return errHandle::line(filename, line_number); }
        SequenceBlock& block = blocks[static_cast<size_t>(current)];
        if (key == "monomers") {
            if (!block.monomers.empty()) block.monomers += ' ';
            block.monomers += value;
        } else if (key == "charges") {
            if (!block.charges.empty()) block.charges += ' ';
            block.charges += value;
        } else if (key == "bonds") {
            if (!block.bonds.empty()) block.bonds += ' ';
            block.bonds += value;
        } else {
            std::fclose(fp);
            return errHandle::unknown("sequence key " + key);
        }
    }
    std::fclose(fp);

    species_list.clear();
    species_list.resize(static_cast<size_t>(p.species_count));
    for (int idx = 0; idx < p.species_count; ++idx) {
        const SequenceBlock& block = blocks[static_cast<size_t>(idx)];
        if (!block.present) return errHandle::message("missing [polymer." + std::to_string(idx) + "] sequence block");
        if (block.monomers.empty()) return errHandle::message("missing monomers for polymer." + std::to_string(idx));
        if (block.charges.empty()) return errHandle::message("missing charges for polymer." + std::to_string(idx));
        const param::polymer_species_input input = p.polymer_species[static_cast<size_t>(idx)];
        if (!build_species_from_sequence(p, input, block, species_list[static_cast<size_t>(idx)])) return false;
    }
    return true;
}

void polymerSpecies::update_counts(param& p)
{
    if (p.species_count != static_cast<int>(species_list.size())) p.species_count = static_cast<int>(species_list.size());
    if (p.polymer_species.size() < species_list.size()) p.polymer_species.resize(species_list.size());

    for (size_t idx = 0; idx < species_list.size(); ++idx) {
        Species& sp = species_list[idx];
        sp.N = static_cast<int>(sp.segments.size());
        if (sp.type_counts.size() != static_cast<size_t>(p.monomer_count())) sp.type_counts.assign(static_cast<size_t>(p.monomer_count()), 0);
        param::polymer_species_input& input = p.polymer_species[idx];
        sp.name = input.name;
        sp.num_chains = input.num_chains;
        input.N = sp.N;
        input.type_counts = sp.type_counts;
    }
    p.recompute_polymer_derived();
}

void polymerSpecies::assign_bond_offsets()
{
    total_bonds = 0;
    for (Species& sp : species_list) {
        sp.bond_offset = total_bonds;
        total_bonds += static_cast<size_t>(std::max(0, sp.N - 1));
    }
}

bool polymerSpecies::same_segment(const Segment& a, const Segment& b) const
{
    const double scale = std::max(1.0, std::max(std::fabs(a.charge), std::fabs(b.charge)));
    return a.type == b.type &&
           std::fabs(a.charge - b.charge) <= 1.0e-12 * scale &&
           std::fabs(a.bond_length - b.bond_length) <= 1.0e-12 * std::max(1.0, std::max(std::fabs(a.bond_length), std::fabs(b.bond_length)));
}

void polymerSpecies::print_sequence(const param& p, const Species& sp) const
{
    std::printf("sequence summary for %s =\n", sp.name.c_str());
    std::printf("%10s %10s %10s %16s %16s %16s\n", "first", "last", "count", "type", "charge", "bond_length");
    int first = 0;
    while (first < static_cast<int>(sp.segments.size())) {
        int last = first;
        while (last + 1 < static_cast<int>(sp.segments.size()) && same_segment(sp.segments[static_cast<size_t>(first)], sp.segments[static_cast<size_t>(last + 1)])) ++last;
        const Segment& segment = sp.segments[static_cast<size_t>(first)];
        const std::string type_name = (segment.type >= 0 && segment.type < static_cast<int>(p.monomer_names.size())) ? p.monomer_names[static_cast<size_t>(segment.type)] : std::to_string(segment.type);
        std::printf("%10d %10d %10d %16s % 16.4e % 16.4e\n", first, last, last - first + 1, type_name.c_str(), segment.charge, segment.bond_length);
        first = last + 1;
    }
}

void polymerSpecies::print(const param& p) const
{
    for (size_t idx = 0; idx < species_list.size(); ++idx) {
        const Species& sp = species_list[idx];
        std::printf("polymer species %zu: name=%s n=%lld N=%d\n", idx, sp.name.c_str(), sp.num_chains, sp.N);
        for (int alpha = 0; alpha < p.monomer_count() && alpha < static_cast<int>(sp.type_counts.size()); ++alpha) {
            std::printf("  count_%s = %d\n", p.monomer_names[static_cast<size_t>(alpha)].c_str(), sp.type_counts[static_cast<size_t>(alpha)]);
        }
        print_sequence(p, sp);
        std::printf("charged segments = %d\n", sp.charged_segment_count());
        std::printf("sum charged z^2 = % 13.4e\n", sp.charge_square_sum());
        std::printf("net chain charge = % 13.4e\n", sp.net_chain_charge());
    }
    std::printf("total polymer charge = % 13.4e\n", total_polymer_charge());
}

void polymerSpecies::build_smoothed_sources(cufftHandle plan, deviceState& d, const param& p, const Complex* omega, const Complex* psi) const
{
    const int size = p.size();
    for (int a = 0; a < p.omega_count(); ++a) {
        copy_and_smooth(plan,
                        d.omega_s + static_cast<size_t>(a) * size,
                        omega + static_cast<size_t>(a) * size,
                        d.Gamma,
                        size);
    }

    for (int alpha = 0; alpha < p.monomer_count(); ++alpha) {
        Complex* W_alpha = d.W_type_s + static_cast<size_t>(alpha) * size;
        clear_complex(W_alpha, size);
        for (int a = 0; a < p.omega_count(); ++a) {
            const Complex coeff = complex_scale(p.omega_gamma(a), p.omega_O(a, alpha));
            add_scaled_complex_const(W_alpha, d.omega_s + static_cast<size_t>(a) * size, coeff, size);
        }
    }

    if (p.electrostatics_enabled()) copy_and_smooth(plan, d.psi_s, psi, d.Gamma, size);
    else clear_complex(d.psi_s, size);
}

void polymerSpecies::build_segment_fields(deviceState& d, const param& p, const Species& sp) const
{
    const int size = p.size();
    for (int s = 0; s < sp.N; ++s) {
        Complex* W = d.W + static_cast<size_t>(s) * size;
        const int alpha = sp.monomer_type(s);
        copy_complex(W, d.W_type_s + static_cast<size_t>(alpha) * size, size);
        const double segment_charge = p.electrostatics_enabled() ? sp.charge(s) : 0.0;
        if (segment_charge != 0.0) add_scaled_i_complex(W, d.psi_s, segment_charge, size);
    }
}

void polymerSpecies::build_exponentials(deviceState& d, const Species& sp, int size) const
{
    for (int s = 0; s < sp.N; ++s) {
        Complex* W = d.W + static_cast<size_t>(s) * size;
        exp_minus_complex(d.exp_mW + static_cast<size_t>(s) * size, W, size);
        exp_complex(d.exp_W + static_cast<size_t>(s) * size, W, size);
    }
}

const Complex* polymerSpecies::exp_minus_for_segment(const deviceState& d, int size, int s) const
{
    return d.exp_mW + static_cast<size_t>(s) * size;
}

const Complex* polymerSpecies::exp_plus_for_segment(const deviceState& d, int size, int s) const
{
    return d.exp_W + static_cast<size_t>(s) * size;
}

void polymerSpecies::propagate_one(cufftHandle plan, Complex* src, Complex* dst, const double* PHI, const Complex* exp_mW, int size) const
{
    fft_apply_real_kernel(plan, src, dst, PHI, size);
    multiply_complex(dst, exp_mW, size);
}

void polymerSpecies::propagate_chain(cufftHandle plan, deviceState& d, const Species& sp, int size) const
{
    copy_complex(d.qF, exp_minus_for_segment(d, size, 0), size);
    copy_complex(d.qB + static_cast<size_t>(sp.N - 1) * size, exp_minus_for_segment(d, size, sp.N - 1), size);
    for (int s = 0; s < sp.N - 1; ++s) {
        const int next = s + 1;
        const double* phi = d.PHI + (sp.bond_offset + static_cast<size_t>(s)) * static_cast<size_t>(size);
        propagate_one(plan, d.qF + static_cast<size_t>(s) * size, d.qF + static_cast<size_t>(next) * size, phi, exp_minus_for_segment(d, size, next), size);
    }
    for (int s = sp.N - 1; s > 0; --s) {
        const int prev = s - 1;
        const double* phi = d.PHI + (sp.bond_offset + static_cast<size_t>(prev)) * static_cast<size_t>(size);
        propagate_one(plan, d.qB + static_cast<size_t>(s) * size, d.qB + static_cast<size_t>(prev) * size, phi, exp_minus_for_segment(d, size, prev), size);
    }
}

Complex polymerSpecies::solve_chain_species_prepared(cufftHandle plan, deviceState& d, const param& p, int species_index) const
{
    const Species& sp = species(species_index);
    const int size = p.size();
    build_segment_fields(d, p, sp);
    build_exponentials(d, sp, size);
    propagate_chain(plan, d, sp, size);
    return mean_complex(d.qF + static_cast<size_t>(sp.N - 1) * size, size);
}

Complex polymerSpecies::solve_chain_species(cufftHandle plan, deviceState& d, const param& p, int species_index, const Complex* omega, const Complex* psi) const
{
    build_smoothed_sources(plan, d, p, omega, psi);
    return solve_chain_species_prepared(plan, d, p, species_index);
}

void polymerSpecies::accumulate_species_densities(deviceState& d, const param& p, const Species& sp, Complex Q) const
{
    const int size = p.size();
    const double chain_density = static_cast<double>(sp.num_chains) / p.volume();
    for (int s = 0; s < sp.N; ++s) {
        const Complex* exp_W = exp_plus_for_segment(d, size, s);
        const int alpha = sp.monomer_type(s);
        Complex* rho_type = d.rho_type + static_cast<size_t>(alpha) * size;
        accumulate_polymer_density_by_chain_density(rho_type, exp_W, d.qF, d.qB, chain_density, Q, 1.0, s, size);
        if (p.electrostatics_enabled()) {
            accumulate_polymer_density_by_chain_density(d.rho_c_poly, exp_W, d.qF, d.qB, chain_density, Q, sp.charge(s), s, size);
        }
    }
}

std::vector<Complex> polymerSpecies::solve_and_compute_densities(cufftHandle plan, deviceState& d, const param& p, const Complex* omega, const Complex* psi) const
{
    const int size = p.size();
    for (int alpha = 0; alpha < p.monomer_count(); ++alpha) clear_complex(d.rho_type + static_cast<size_t>(alpha) * size, size);
    clear_complex(d.rho_c_poly, size);
    clear_complex(d.rho_c_ion, size);
    clear_complex(d.rho_c, size);

    build_smoothed_sources(plan, d, p, omega, psi);

    std::vector<Complex> Qs;
    Qs.reserve(species_list.size());
    for (int idx = 0; idx < static_cast<int>(species_list.size()); ++idx) {
        const Complex Q = solve_chain_species_prepared(plan, d, p, idx);
        Qs.push_back(Q);
        accumulate_species_densities(d, p, species_list[static_cast<size_t>(idx)], Q);
    }

    clear_complex(d.rho_p, size);
    for (int alpha = 0; alpha < p.monomer_count(); ++alpha) {
        Complex* rho_alpha = d.rho_type + static_cast<size_t>(alpha) * size;
        fft_smooth(plan, rho_alpha, d.Gamma, size);
        add_complex(d.rho_p, rho_alpha, size);
    }
    if (p.electrostatics_enabled()) {
        fft_smooth(plan, d.rho_c_poly, d.Gamma, size);
        add_complex(d.rho_c, d.rho_c_poly, size);
    }
    return Qs;
}
