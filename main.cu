#include "complexLangevin.h"
#include "fieldIO.h"
#include "spectral.h"

#include <algorithm>
#include <chrono>
#include <cstdio>

int main(int argc, char** argv)
{
    const auto start_time = std::chrono::steady_clock::now();
    parse::options options;
    if (!parse::args(argc, argv, options)) return EXIT_FAILURE;

    param p;
    if (!p.read(options.param_file)) return EXIT_FAILURE;
    if (options.sequence_override) p.sequence_file = options.sequence_file;

    fieldIO io;
    if (!io.prepare(options, p)) return EXIT_FAILURE;

    polymerSpecies poly;
    if (!poly.setup(p)) return EXIT_FAILURE;
    if (!p.validate()) return EXIT_FAILURE;

    ion ions;
    if (!ions.build(p, poly.total_polymer_charge())) return EXIT_FAILURE;

    if (options.screenprint) {
        std::printf("parameter file = %s\n", options.param_file);
        if (!io.run_path().empty()) std::printf("run directory = %s\n", io.run_path().c_str());
        std::printf("omega input file = %s\n", io.input_omega_file().c_str());
        if (p.electrostatics_enabled()) std::printf("psi input file = %s\n", io.input_psi_file().c_str());
        else std::printf("psi input file = ignored (electrostatics disabled)\n");
        std::printf("omega output file = %s\n", p.omega_file.c_str());
        std::printf("psi output file = %s\n", p.psi_file.c_str());
        if (p.Ntmp > 0) {
            std::printf("tmp omega checkpoint file = %s\n", p.tmp_omega_file.c_str());
            std::printf("tmp psi checkpoint file = %s\n", p.tmp_psi_file.c_str());
        }
        std::printf("observable file = %s\n", p.observable_file.c_str());
        std::printf("sequence file = %s\n", p.sequence_file.c_str());
        poly.print(p);
        p.print();
        ions.print(p);
    }

    const int size = p.size();
    std::vector<Complex> h_omega(static_cast<size_t>(p.omega_count()) * size);
    std::vector<Complex> h_psi(size);
    std::vector<Complex> h_rho_type(static_cast<size_t>(p.monomer_count()) * size);
    std::vector<Complex> h_rho_p(size);
    std::vector<Complex> h_rho_c(size);
    std::vector<Complex> h_rho_c_poly(size);
    std::vector<Complex> h_rho_c_ion(size);
    std::vector<Complex> h_ion_rho(static_cast<size_t>(size) * ions.count());
    std::vector<double> h_Gamma;
    std::vector<double> h_PHI;
    std::vector<double> h_k2;
    std::vector<double> h_ion_kernels;

    spectral::build_kernels(p, poly, h_Gamma, h_PHI, h_k2);
    ions.build_kernels(p, h_k2, h_ion_kernels);

    if (!io.read_omega_fields(p, h_omega)) return EXIT_FAILURE;
    if (p.electrostatics_enabled()) {
        if (!io.read_psi_field(p, h_psi)) return EXIT_FAILURE;
    } else {
        std::fill(h_psi.begin(), h_psi.end(), complex_zero());
    }
    if (!finite_complex_vectors({{"input omega", &h_omega}, {"input psi", &h_psi}})) return EXIT_FAILURE;

    complexLangevin cl(p, poly, ions);
    if (!cl.initialize(h_omega, h_psi, h_Gamma, h_PHI, h_k2, h_ion_kernels)) return EXIT_FAILURE;
    if (!cl.run(options.screenprint)) return EXIT_FAILURE;
    if (!cl.final_densities()) return EXIT_FAILURE;
    cl.copy_results(h_omega, h_psi, h_rho_type, h_rho_p, h_rho_c, h_rho_c_poly, h_rho_c_ion, h_ion_rho);

    if (!finite_complex_vectors({{"output omega", &h_omega},
                                {"output psi", &h_psi},
                                {"output rho_type", &h_rho_type},
                                {"output rho_p", &h_rho_p},
                                {"output rho_c", &h_rho_c},
                                {"output rho_c_poly", &h_rho_c_poly},
                                {"output rho_c_ion", &h_rho_c_ion},
                                {"output ion_rho", &h_ion_rho}})) return EXIT_FAILURE;

    const bool success = io.write_all(p, ions, h_omega, h_psi, h_rho_type, h_rho_p, h_rho_c, h_rho_c_poly, h_rho_c_ion, h_ion_rho);
    const auto end_time = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(end_time - start_time).count();
    std::printf("total time = % .4e s\n", seconds);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
