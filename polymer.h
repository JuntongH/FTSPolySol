#pragma once

#include "param.h"
#include "deviceState.h"

class polymerSpecies {
public:
    struct Species {
        std::string name;
        long long num_chains = 0;
        int N = 0;
        std::vector<int> type_counts;
        size_t bond_offset = 0;
        std::vector<Segment> segments;

        int monomer_type(int s) const { return segments[static_cast<size_t>(s)].type; }
        double charge(int s) const { return segments[static_cast<size_t>(s)].charge; }
        double bond_length(int s) const { return segments[static_cast<size_t>(s)].bond_length; }
        int charged_segment_count() const;
        double charge_square_sum() const;
        double net_chain_charge() const;
    };

    bool setup(param& p);
    int count() const;
    const Species& species(int idx) const;
    const std::vector<Species>& all() const;
    long long total_chains() const;
    size_t bond_kernel_count() const;
    double total_polymer_charge() const;
    void print(const param& p) const;

    void build_smoothed_sources(cufftHandle plan,
                                deviceState& d,
                                const param& p,
                                const Complex* omega,
                                const Complex* psi) const;

    Complex solve_chain_species(cufftHandle plan,
                                deviceState& d,
                                const param& p,
                                int species_index,
                                const Complex* omega,
                                const Complex* psi) const;

    Complex solve_chain_species_prepared(cufftHandle plan,
                                         deviceState& d,
                                         const param& p,
                                         int species_index) const;

    std::vector<Complex> solve_and_compute_densities(cufftHandle plan,
                                                     deviceState& d,
                                                     const param& p,
                                                     const Complex* omega,
                                                     const Complex* psi) const;

private:
    std::vector<Species> species_list;
    size_t total_bonds = 0;

    bool read_block_sequence_file(const char* filename, param& p);
    void update_counts(param& p);
    void assign_bond_offsets();
    bool same_segment(const Segment& a, const Segment& b) const;
    void print_sequence(const param& p, const Species& sp) const;
    void build_segment_fields(deviceState& d, const param& p, const Species& sp) const;
    void build_exponentials(deviceState& d, const Species& sp, int size) const;
    void propagate_chain(cufftHandle plan, deviceState& d, const Species& sp, int size) const;
    void propagate_one(cufftHandle plan, Complex* src, Complex* dst, const double* PHI, const Complex* exp_mW, int size) const;
    void accumulate_species_densities(deviceState& d, const param& p, const Species& sp, Complex Q) const;
    const Complex* exp_minus_for_segment(const deviceState& d, int size, int s) const;
    const Complex* exp_plus_for_segment(const deviceState& d, int size, int s) const;
};

