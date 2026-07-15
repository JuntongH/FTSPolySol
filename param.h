#pragma once

#include "parse.h"
#include "interactionMode.h"

class param {
public:
    struct polymer_species_input {
        std::string name;
        long long num_chains = 0;
        int N = 0;
        std::vector<int> type_counts;
        bool section_seen = false;
        bool name_set = false;
        bool num_chains_set = false;
        bool N_set = false;
    };

    struct counterion_input {
        bool enabled = false;
        bool auto_number = true;
        long long number = 0;
        double valence = 0.0;
        double smearing_length = UnsetSmearing;
        std::string name = "counterion";
    };

    struct salt_input {
        bool enabled = false;
        long long num_salt = 0;
        double cation_valence = 0.0;
        double anion_valence = 0.0;
        int cation_stoich = 1;
        int anion_stoich = 1;
        double smearing_length = UnsetSmearing;
        bool combine_identical_species = true;
        std::string cation_name = "salt_cation";
        std::string anion_name = "salt_anion";
    };

    struct stress_input {
        bool enabled = false;
        bool component_xx = false;
        bool component_yy = false;
        bool component_zz = false;
        bool component_xy = false;
        bool component_xz = false;
        bool component_yz = false;
        bool isotropic_pressure = false;

        bool any_component() const
        {
            return component_xx || component_yy || component_zz ||
                   component_xy || component_xz || component_yz;
        }

        bool any_observable() const
        {
            return any_component() || isotropic_pressure;
        }
    };

    int species_count = 1;
    int Nmax = 0;
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    int Nt = 0;
    int Nss = 0;
    int Ns = 0;
    int Ntmp = 0;
    int Nso = -1;
    int N_stat = -1;
    int max_nan_restarts = 20;
    long long total_polymer_chains = 0;
    long long total_segments = 0;
    std::vector<long long> total_type_segments;

    int monomer_type_count = 2;
    bool monomer_count_set = false;
    bool monomer_names_set = false;
    bool chi_matrix_set = false;
    std::vector<std::string> monomer_names;
    std::vector<double> chi_matrix;
    double eigen_tol_rel = 1.0e-12;
    interactionMode modes;

    double a = 0.0;
    double dL = 0.0;
    double xi = 0.0;
    double lB = 0.0;
    double L1 = 0.0;
    double L2 = 0.0;
    double L3 = 0.0;
    double rho0 = 0.0;
    double dt = 0.0;
    double lambdaPsi = 0.0;
    double dV = 0.0;
    double default_smearing_length = UnsetSmearing;
    bool noise = true;
    bool ions_enabled = false;
    std::vector<double> lambda_omega_input;
    std::vector<double> lambda_omega;

    std::string sequence_file;
    std::string omega_file = "omega.rf";             // active omega-mode fields
    std::string psi_file = "psi.rf";
    std::string tmp_omega_file = "tmp_omega.rf";
    std::string tmp_psi_file = "tmp_psi.rf";
    std::string rho_type_file = "rho_type.rf";   // one field per monomer type
    std::string rhop_file = "rhop.rf";
    std::string rhoc_file = "rhoc.rf";
    std::string rhoc_poly_file = "rhoc_poly.rf";
    std::string rhoc_ion_file = "rhoc_ion.rf";
    std::string ion_density_prefix = "rho_ion";
    std::string step_omega_file;
    std::string step_psi_file;
    std::string step_rho_type_file;
    std::string step_rhop_file;
    std::string step_rhoc_file;
    std::string step_rhoc_poly_file;
    std::string step_rhoc_ion_file;
    std::string step_ion_density_prefix;
    std::string observable_file = "observables.dat";
    std::string drift_max_file = "drift_max.dat";

    std::vector<polymer_species_input> polymer_species;
    bool species_count_set = false;
    counterion_input counterions;
    salt_input salt;
    stress_input stress;

    bool read(const char* filename);
    bool validate() const;
    bool validate_modes() const;
    void print() const;
    int size() const;
    double volume() const;
    bool electrostatics_enabled() const;
    int index(int x, int y, int z) const;
    int monomer_index(const std::string& name) const;
    int omega_count() const;
    int monomer_count() const;
    double omega_coefficient(int active_index) const;
    double omega_mobility(int active_index) const;
    Complex omega_gamma(int active_index) const;
    double omega_O(int active_index, int alpha) const;
    double short_range_self_pressure() const;
    bool legacy_ab_field_order() const;
    void recompute_polymer_derived();

private:
    bool apply(const std::string& section, const std::string& key, const std::string& value);
    bool unknown(const std::string& section, const std::string& key) const;
    bool derive();
    bool finalize_interactions();
    bool finalize_mobilities();
    polymer_species_input& ensure_polymer_species(int index);
};
