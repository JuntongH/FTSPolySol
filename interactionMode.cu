#include "interactionMode.h"
#include "errHandle.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

bool interactionMode::build(int monomer_type_count,
                            const std::vector<std::string>& monomer_names,
                            double xi,
                            const std::vector<double>& chi_matrix,
                            double eigen_tol_rel)
{
    if (!validate_input(monomer_type_count, monomer_names, xi, chi_matrix, eigen_tol_rel)) return false;

    M_ = monomer_type_count;
    monomer_names_ = monomer_names;
    xi_ = xi;
    chi_ = chi_matrix;
    eigen_tol_rel_ = eigen_tol_rel;
    eigen_tol_abs_ = 0.0;
    reconstruction_error_ = 0.0;
    orthogonality_error_ = 0.0;
    K_.clear();
    eigenvalues_.clear();
    O_.clear();
    active_modes_.clear();

    const std::size_t mm = static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_);

    double max_chi = 0.0;
    for (double x : chi_) max_chi = std::max(max_chi, std::abs(x));
    const double sym_tol = 1.0e-12 * std::max(1.0, max_chi);

    for (int alpha = 0; alpha < M_; ++alpha) {
        const std::size_t aa = static_cast<std::size_t>(alpha) * M_ + alpha;
        if (std::abs(chi_[aa]) > sym_tol) return errHandle::invalid("chi diagonal entries must be zero");
        chi_[aa] = 0.0;
        for (int beta = alpha + 1; beta < M_; ++beta) {
            const std::size_t ab = static_cast<std::size_t>(alpha) * M_ + beta;
            const std::size_t ba = static_cast<std::size_t>(beta) * M_ + alpha;
            if (std::abs(chi_[ab] - chi_[ba]) > sym_tol) return errHandle::invalid("chi matrix must be symmetric");
            const double avg = 0.5 * (chi_[ab] + chi_[ba]);
            chi_[ab] = avg;
            chi_[ba] = avg;
        }
    }

    if (!build_K()) return false;
    if (K_.size() != mm) return errHandle::invalid("K matrix size");
    if (!diagonalize()) return false;
    if (!compute_checks()) return false;

    if (reconstruction_error_ > 1.0e-10) return errHandle::message("interaction eigenmode reconstruction check failed");
    if (orthogonality_error_ > 1.0e-10) return errHandle::message("interaction eigenvector orthogonality check failed");

    return true;
}

bool interactionMode::validate_input(int monomer_type_count,
                                     const std::vector<std::string>& monomer_names,
                                     double xi,
                                     const std::vector<double>& chi_matrix,
                                     double eigen_tol_rel) const
{
    if (monomer_type_count <= 0) return errHandle::invalid("[monomers] M");
    if (monomer_names.size() != static_cast<std::size_t>(monomer_type_count)) {
        return errHandle::invalid("[monomers] names length must equal M");
    }
    for (int i = 0; i < monomer_type_count; ++i) {
        if (monomer_names[static_cast<std::size_t>(i)].empty()) return errHandle::invalid("empty monomer name");
        for (int j = i + 1; j < monomer_type_count; ++j) {
            if (monomer_names[static_cast<std::size_t>(i)] == monomer_names[static_cast<std::size_t>(j)]) {
                return errHandle::invalid("duplicate monomer name " + monomer_names[static_cast<std::size_t>(i)]);
            }
        }
    }
    if (!std::isfinite(xi) || xi < 0.0) return errHandle::invalid("xi");
    if (!std::isfinite(eigen_tol_rel) || eigen_tol_rel < 0.0) return errHandle::invalid("eigen_tol_rel");
    const std::size_t mm = static_cast<std::size_t>(monomer_type_count) * static_cast<std::size_t>(monomer_type_count);
    if (chi_matrix.size() != mm) return errHandle::invalid("[interactions] chi matrix size must be M*M");
    for (double x : chi_matrix) {
        if (!std::isfinite(x)) return errHandle::invalid("chi matrix entry");
    }
    return true;
}

bool interactionMode::build_K()
{
    const std::size_t mm = static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_);
    K_.assign(mm, 0.0);
    for (int alpha = 0; alpha < M_; ++alpha) {
        for (int beta = 0; beta < M_; ++beta) {
            const std::size_t ab = static_cast<std::size_t>(alpha) * M_ + beta;
            K_[ab] = (alpha == beta) ? xi_ : xi_ + chi_[ab];
        }
    }
    return true;
}

double interactionMode::matrix_frobenius(const std::vector<double>& a)
{
    long double sum = 0.0L;
    for (double x : a) sum += static_cast<long double>(x) * static_cast<long double>(x);
    return std::sqrt(static_cast<double>(sum));
}

double interactionMode::offdiag_frobenius(const std::vector<double>& a, int n)
{
    long double sum = 0.0L;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double x = a[static_cast<std::size_t>(i) * n + j];
            sum += static_cast<long double>(x) * static_cast<long double>(x);
        }
    }
    return std::sqrt(static_cast<double>(sum));
}

bool interactionMode::jacobi_symmetric_eigen(const std::vector<double>& input,
                                             int n,
                                             std::vector<double>& eigenvalues,
                                             std::vector<double>& eigenvectors)
{
    if (n <= 0) return errHandle::invalid("interaction matrix dimension");
    if (input.size() != static_cast<std::size_t>(n) * static_cast<std::size_t>(n)) {
        return errHandle::invalid("interaction matrix size");
    }

    std::vector<double> a = input;
    eigenvectors.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) eigenvectors[static_cast<std::size_t>(i) * n + i] = 1.0;

    if (n == 1) {
        eigenvalues.assign(1, a[0]);
        return true;
    }

    const double scale = std::max(1.0, matrix_frobenius(input));
    const double tol = 1.0e-14 * scale;
    const int max_sweeps = std::max(50, 100 * n * n);
    bool converged = false;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        if (offdiag_frobenius(a, n) <= tol) {
            converged = true;
            break;
        }
        bool rotated = false;
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const std::size_t pp = static_cast<std::size_t>(p) * n + p;
                const std::size_t qq = static_cast<std::size_t>(q) * n + q;
                const std::size_t pq = static_cast<std::size_t>(p) * n + q;
                const double apq = a[pq];
                if (std::abs(apq) <= std::numeric_limits<double>::epsilon() * scale) continue;

                const double app = a[pp];
                const double aqq = a[qq];
                const double tau = (aqq - app) / (2.0 * apq);
                const double t = (tau >= 0.0)
                                     ? 1.0 / (tau + std::sqrt(1.0 + tau * tau))
                                     : -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = t * c;

                for (int k = 0; k < n; ++k) {
                    if (k == p || k == q) continue;
                    const std::size_t kp = static_cast<std::size_t>(k) * n + p;
                    const std::size_t kq = static_cast<std::size_t>(k) * n + q;
                    const double akp = a[kp];
                    const double akq = a[kq];
                    const double new_kp = c * akp - s * akq;
                    const double new_kq = s * akp + c * akq;
                    a[kp] = new_kp;
                    a[static_cast<std::size_t>(p) * n + k] = new_kp;
                    a[kq] = new_kq;
                    a[static_cast<std::size_t>(q) * n + k] = new_kq;
                }

                a[pp] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                a[qq] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[pq] = 0.0;
                a[static_cast<std::size_t>(q) * n + p] = 0.0;

                for (int k = 0; k < n; ++k) {
                    const std::size_t kp = static_cast<std::size_t>(k) * n + p;
                    const std::size_t kq = static_cast<std::size_t>(k) * n + q;
                    const double vkp = eigenvectors[kp];
                    const double vkq = eigenvectors[kq];
                    eigenvectors[kp] = c * vkp - s * vkq;
                    eigenvectors[kq] = s * vkp + c * vkq;
                }
                rotated = true;
            }
        }
        if (!rotated) {
            converged = true;
            break;
        }
    }

    if (!converged && offdiag_frobenius(a, n) > 1.0e-10 * scale) {
        return errHandle::message("Jacobi diagonalization of interaction matrix did not converge");
    }

    eigenvalues.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) eigenvalues[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i) * n + i];
    return true;
}

bool interactionMode::diagonalize()
{
    std::vector<double> raw_eigenvalues;
    std::vector<double> raw_eigenvectors; // V_{alpha,l}, row-major by alpha, column by l.
    if (!jacobi_symmetric_eigen(K_, M_, raw_eigenvalues, raw_eigenvectors)) return false;

    std::vector<int> order(static_cast<std::size_t>(M_));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const double da = raw_eigenvalues[static_cast<std::size_t>(a)];
        const double db = raw_eigenvalues[static_cast<std::size_t>(b)];
        if (da != db) return da > db;
        return a < b;
    });

    eigenvalues_.assign(static_cast<std::size_t>(M_), 0.0);
    O_.assign(static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_), 0.0);

    for (int ell = 0; ell < M_; ++ell) {
        const int raw = order[static_cast<std::size_t>(ell)];
        eigenvalues_[static_cast<std::size_t>(ell)] = raw_eigenvalues[static_cast<std::size_t>(raw)];
        for (int alpha = 0; alpha < M_; ++alpha) {
            O_[static_cast<std::size_t>(ell) * M_ + alpha] =
                raw_eigenvectors[static_cast<std::size_t>(alpha) * M_ + raw];
        }

        int pivot = 0;
        double pivot_abs = std::abs(O_[static_cast<std::size_t>(ell) * M_]);
        for (int alpha = 1; alpha < M_; ++alpha) {
            const double mag = std::abs(O_[static_cast<std::size_t>(ell) * M_ + alpha]);
            if (mag > pivot_abs) {
                pivot_abs = mag;
                pivot = alpha;
            }
        }
        if (O_[static_cast<std::size_t>(ell) * M_ + pivot] < 0.0) {
            for (int alpha = 0; alpha < M_; ++alpha) O_[static_cast<std::size_t>(ell) * M_ + alpha] *= -1.0;
        }
    }

    double max_abs_d = 0.0;
    for (double d : eigenvalues_) max_abs_d = std::max(max_abs_d, std::abs(d));
    eigen_tol_abs_ = eigen_tol_rel_ * std::max(1.0, max_abs_d);

    active_modes_.clear();
    for (int ell = 0; ell < M_; ++ell) {
        if (std::abs(eigenvalues_[static_cast<std::size_t>(ell)]) > eigen_tol_abs_) active_modes_.push_back(ell);
    }
    return true;
}

bool interactionMode::compute_checks()
{
    const std::size_t mm = static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_);
    std::vector<double> reconstructed(mm, 0.0);
    for (int alpha = 0; alpha < M_; ++alpha) {
        for (int beta = 0; beta < M_; ++beta) {
            double value = 0.0;
            for (int ell = 0; ell < M_; ++ell) {
                const double oa = O_[static_cast<std::size_t>(ell) * M_ + alpha];
                const double ob = O_[static_cast<std::size_t>(ell) * M_ + beta];
                value += oa * eigenvalues_[static_cast<std::size_t>(ell)] * ob;
            }
            reconstructed[static_cast<std::size_t>(alpha) * M_ + beta] = value;
        }
    }

    std::vector<double> residual(mm, 0.0);
    for (std::size_t i = 0; i < mm; ++i) residual[i] = K_[i] - reconstructed[i];
    reconstruction_error_ = matrix_frobenius(residual) / std::max(1.0, matrix_frobenius(K_));

    std::vector<double> orth(mm, 0.0);
    for (int ell = 0; ell < M_; ++ell) {
        for (int m = 0; m < M_; ++m) {
            double dot = 0.0;
            for (int alpha = 0; alpha < M_; ++alpha) {
                dot += O_[static_cast<std::size_t>(ell) * M_ + alpha] *
                       O_[static_cast<std::size_t>(m) * M_ + alpha];
            }
            orth[static_cast<std::size_t>(ell) * M_ + m] = dot - (ell == m ? 1.0 : 0.0);
        }
    }
    orthogonality_error_ = matrix_frobenius(orth);
    return std::isfinite(reconstruction_error_) && std::isfinite(orthogonality_error_);
}

bool interactionMode::valid() const
{
    const std::size_t mm = static_cast<std::size_t>(std::max(0, M_)) * static_cast<std::size_t>(std::max(0, M_));
    return M_ > 0 &&
           monomer_names_.size() == static_cast<std::size_t>(M_) &&
           chi_.size() == mm &&
           K_.size() == mm &&
           eigenvalues_.size() == static_cast<std::size_t>(M_) &&
           O_.size() == mm &&
           std::isfinite(eigen_tol_rel_) && eigen_tol_rel_ >= 0.0 &&
           std::isfinite(reconstruction_error_) && reconstruction_error_ <= 1.0e-10 &&
           std::isfinite(orthogonality_error_) && orthogonality_error_ <= 1.0e-10;
}

int interactionMode::monomer_index(const std::string& name) const
{
    for (int i = 0; i < static_cast<int>(monomer_names_.size()); ++i) {
        if (monomer_names_[static_cast<std::size_t>(i)] == name) return i;
    }
    return -1;
}

double interactionMode::O(int mode, int alpha) const
{
    return O_.at(static_cast<std::size_t>(mode) * M_ + alpha);
}


Complex interactionMode::active_gamma(int active_index) const
{
    const double d = active_eigenvalue(active_index);
    return d > 0.0 ? make_complex(0.0, 1.0) : make_complex(1.0, 0.0);
}

const char* interactionMode::gamma_symbol(int mode) const
{
    if (mode < 0 || mode >= M_) return "-";
    const bool active = std::find(active_modes_.begin(), active_modes_.end(), mode) != active_modes_.end();
    if (!active) return "-";
    return eigenvalues_[static_cast<std::size_t>(mode)] > 0.0 ? "i" : "1";
}

void interactionMode::print_mode_table(FILE* out) const
{
    std::fprintf(out, "monomer types (M = %d):", M_);
    for (const std::string& name : monomer_names_) std::fprintf(out, " %s", name.c_str());
    std::fprintf(out, "\n");

    std::fprintf(out, "xi = % .17g\n", xi_);
    std::fprintf(out, "eigen_tol_rel = %.17g\n", eigen_tol_rel_);
    std::fprintf(out, "eigen_tol_abs = %.17g\n", eigen_tol_abs_);

    if (chi_.size() == static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_)) {
        std::fprintf(out, "chi matrix:\n");
        for (int alpha = 0; alpha < M_; ++alpha) {
            std::fprintf(out, "  ");
            for (int beta = 0; beta < M_; ++beta) {
                std::fprintf(out, " % .10e", chi_[static_cast<std::size_t>(alpha) * M_ + beta]);
            }
            std::fprintf(out, "\n");
        }
    }

    if (K_.size() == static_cast<std::size_t>(M_) * static_cast<std::size_t>(M_)) {
        std::fprintf(out, "K matrix, K_aa = xi and K_ab = xi + chi[a,b]:\n");
        for (int alpha = 0; alpha < M_; ++alpha) {
            std::fprintf(out, "  ");
            for (int beta = 0; beta < M_; ++beta) {
                std::fprintf(out, " % .10e", K_[static_cast<std::size_t>(alpha) * M_ + beta]);
            }
            std::fprintf(out, "\n");
        }
    }

    std::fprintf(out, "interaction eigenmodes, K = O^T D O:\n");
    std::fprintf(out, "%-6s %-8s %-18s %-7s", "mode", "active", "d", "gamma");
    for (const std::string& name : monomer_names_) std::fprintf(out, " %16s", name.c_str());
    std::fprintf(out, "\n");

    for (int ell = 0; ell < M_; ++ell) {
        const double d = eigenvalues_[static_cast<std::size_t>(ell)];
        const bool active = std::find(active_modes_.begin(), active_modes_.end(), ell) != active_modes_.end();
        std::fprintf(out, "%-6d %-8s % .10e %-7s", ell, active ? "yes" : "no", d, gamma_symbol(ell));
        for (int alpha = 0; alpha < M_; ++alpha) {
            std::fprintf(out, " % 16.10e", O_[static_cast<std::size_t>(ell) * M_ + alpha]);
        }
        std::fprintf(out, "\n");
    }

    std::fprintf(out, "active mode set I = {");
    for (std::size_t i = 0; i < active_modes_.size(); ++i) {
        if (i > 0) std::fprintf(out, ", ");
        std::fprintf(out, "%d", active_modes_[i]);
    }
    std::fprintf(out, "}\n");
    std::fprintf(out, "checks: ||K - O^T D O||_F/max(1,||K||_F) = %.3e, ||O O^T - I||_F = %.3e\n",
                 reconstruction_error_, orthogonality_error_);
}
