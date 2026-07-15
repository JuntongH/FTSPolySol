#pragma once

#include "param.h"
#include "ion.h"
#include "trilinear.h"
#include <sys/stat.h>
#include <sys/types.h>

class fieldIO {
public:
    bool prepare(const parse::options& options, param& p);
    const std::string& run_path() const;
    const std::string& input_omega_file() const;
    const std::string& input_psi_file() const;

    bool read_omega_fields(const param& p, std::vector<Complex>& omega) const;
    bool read_psi_field(const param& p, std::vector<Complex>& psi) const;
    bool write_all(const param& p, const ion& ions, const std::vector<Complex>& omega, const std::vector<Complex>& psi, const std::vector<Complex>& rho_type, const std::vector<Complex>& rho_p, const std::vector<Complex>& rho_c, const std::vector<Complex>& rho_c_poly, const std::vector<Complex>& rho_c_ion, const std::vector<Complex>& ion_rho) const;

    static bool read_omega_fields_exact(const param& p, const std::string& filename, std::vector<Complex>& omega);
    static bool read_complex_field_exact(const param& p, const std::string& filename, std::vector<Complex>& field);
    static bool write_omega_fields(const param& p, const std::string& filename, const std::vector<Complex>& omega);
    static bool write_type_fields(const param& p, const std::string& filename, const std::vector<Complex>& fields, int field_count);
    static bool write_complex_field(const param& p, const std::string& filename, const std::vector<Complex>& field);
    static bool step_outputs_enabled(const param& p);
    static bool write_step_outputs(const param& p, const ion& ions, int output_index, const std::vector<Complex>& omega, const std::vector<Complex>& psi, const std::vector<Complex>& rho_type, const std::vector<Complex>& rho_p, const std::vector<Complex>& rho_c, const std::vector<Complex>& rho_c_poly, const std::vector<Complex>& rho_c_ion, const std::vector<Complex>& ion_rho);

private:
    std::string run;
    std::string input_omega;
    std::string input_psi;
    bool screenprint = false;
    bool zero_missing_default_omega = false;
    mutable bool source_mesh_known = false;
    mutable int source_m1 = 0;
    mutable int source_m2 = 0;
    mutable int source_m3 = 0;

    static bool path_exists_as_directory(const std::string& path);
    static bool make_directory(const std::string& path);
    static bool make_directories(const std::string& path);
    static bool copy_file(const std::string& source, const std::string& target);
    static bool copy_case_files(const parse::options& options, const param& p, const std::string& run_path);
    static bool valid_run_id(const std::string& id);
    static std::string join_path(const std::string& a, const std::string& b);
    static std::string directory_path(const std::string& path);
    static std::string basename_path(const std::string& path);
    static std::string sequence_directory_name(const std::string& pattern);
    static std::string run_output_path(const std::string& run_path, const std::string& name);
    static bool prepare_step_output_template(std::string& pattern);
    static bool prepare_step_output_templates(param& p);
    static bool infer_mesh(size_t rows, const param& p, int& m1, int& m2, int& m3);
    static bool write_complex_field(const param& p, const std::string& filename, const Complex* field, size_t size);
    static bool write_ion_densities_with_prefix(const param& p, const ion& ions, const std::string& prefix, const std::vector<Complex>& ion_rho);
    static bool write_ion_densities(const param& p, const ion& ions, const std::vector<Complex>& ion_rho);
};
