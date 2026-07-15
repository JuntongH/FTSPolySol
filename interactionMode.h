#pragma once

#include "common.h"

#include <cstdio>
#include <string>
#include <vector>

class interactionMode {
public:
    bool build(int monomer_type_count,
               const std::vector<std::string>& monomer_names,
               double xi,
               const std::vector<double>& chi_matrix,
               double eigen_tol_rel);

    bool valid() const;
    void print_mode_table(FILE* out = stdout) const;

    int monomer_count() const { return M_; }
    int active_count() const { return static_cast<int>(active_modes_.size()); }
    int active_mode(int i) const { return active_modes_.at(static_cast<std::size_t>(i)); }
    int monomer_index(const std::string& name) const;

    const std::vector<std::string>& monomer_names() const { return monomer_names_; }
    const std::vector<double>& chi_matrix() const { return chi_; }
    const std::vector<double>& K_matrix() const { return K_; }
    const std::vector<double>& eigenvalues() const { return eigenvalues_; }
    const std::vector<double>& eigenvectors() const { return O_; }
    const std::vector<int>& active_modes() const { return active_modes_; }

    double xi() const { return xi_; }
    double eigen_tol_rel() const { return eigen_tol_rel_; }
    double eigen_tol_abs() const { return eigen_tol_abs_; }
    double reconstruction_error() const { return reconstruction_error_; }
    double orthogonality_error() const { return orthogonality_error_; }
    double O(int mode, int alpha) const;
    double active_O(int active_index, int alpha) const { return O(active_mode(active_index), alpha); }
    double eigenvalue(int mode) const { return eigenvalues_.at(static_cast<std::size_t>(mode)); }
    double active_eigenvalue(int active_index) const { return eigenvalue(active_mode(active_index)); }
    double active_inv_abs_eigenvalue(int active_index) const { return 1.0 / std::abs(active_eigenvalue(active_index)); }
    Complex active_gamma(int active_index) const;
    const char* gamma_symbol(int mode) const;

private:
    static double matrix_frobenius(const std::vector<double>& a);
    static double offdiag_frobenius(const std::vector<double>& a, int n);
    static bool jacobi_symmetric_eigen(const std::vector<double>& input,
                                       int n,
                                       std::vector<double>& eigenvalues,
                                       std::vector<double>& eigenvectors);

    bool validate_input(int monomer_type_count,
                        const std::vector<std::string>& monomer_names,
                        double xi,
                        const std::vector<double>& chi_matrix,
                        double eigen_tol_rel) const;
    bool build_K();
    bool diagonalize();
    bool compute_checks();

    int M_ = 0;
    std::vector<std::string> monomer_names_;
    double xi_ = 0.0;
    std::vector<double> chi_;          // chi_{alpha beta}, row-major, M*M.
    std::vector<double> K_;            // K_{alpha beta}, row-major, M*M.
    std::vector<double> eigenvalues_;  // d_l, sorted by decreasing d_l.
    std::vector<double> O_;            // O_{l alpha}, row-major by mode l.
    std::vector<int> active_modes_;    // l with |d_l| > eigen_tol_abs_.
    double eigen_tol_rel_ = 1.0e-12;
    double eigen_tol_abs_ = 0.0;
    double reconstruction_error_ = 0.0;
    double orthogonality_error_ = 0.0;
};
