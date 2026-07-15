#pragma once

#include "polymer.h"
#include "ion.h"
#include "observable.h"

class complexLangevin {
public:
    complexLangevin(param& p, polymerSpecies& poly, ion& ions);
    ~complexLangevin();

    bool initialize(const std::vector<Complex>& omega,
                    const std::vector<Complex>& psi,
                    const std::vector<double>& Gamma,
                    const std::vector<double>& PHI,
                    const std::vector<double>& k2,
                    const std::vector<double>& ion_kernels);
    bool run(bool screenprint);
    bool final_densities();

    void copy_results(std::vector<Complex>& omega,
                      std::vector<Complex>& psi,
                      std::vector<Complex>& rho_type,
                      std::vector<Complex>& rho_p,
                      std::vector<Complex>& rho_c,
                      std::vector<Complex>& rho_c_poly,
                      std::vector<Complex>& rho_c_ion,
                      std::vector<Complex>& ion_rho);

private:
    enum class StepStatus {
        ok,
        failed,
        nonfinite
    };

    param& p;
    polymerSpecies& poly;
    ion& ions;
    observable obs;
    deviceState d;
    cufftHandle plan = 0;
    bool ready = false;

    FILE* observable_output = nullptr;
    FILE* drift_output = nullptr;
    std::vector<double> observable_sum;
    long long observable_count = 0;
    bool checkpoint_valid = false;
    int checkpoint_step = -1;
    int restart_count = 0;

    std::vector<Complex> compute_densities(const Complex* omega, const Complex* psi);
    bool finite_chain_partitions(const std::vector<Complex>& Qp, const char* label, int step) const;
    void compute_drifts(const Complex* omega_term, Complex* drift);
    void generate_noises();
    void compute_psi_drift(const Complex* psi_field, const Complex* rho_source, Complex* drift);
    void remove_psi_zero_mode(Complex* psi_field);
    double corrector_drift_max(Complex* scratch, const Complex* first, const Complex* second, double mobility);
    void predictor_step();
    void corrector_step();

    bool open_observables();
    bool open_drift_statistics();
    void write_observable_header();
    void write_drift_header();
    void write_observable_record(int step, const observableRecord& record);
    bool accumulate_observable_record(int step, const observableRecord& record);
    bool write_observable_average();
    bool drift_statistics_enabled() const;
    bool write_drift_max(int step);
    bool reset_drift_statistics_after_restart(int failed_step);
    bool checkpointing_enabled() const;
    bool write_step_outputs(int step, bool& nonfinite);
    bool write_checkpoint(int step, bool& nonfinite);
    bool read_checkpoint(int failed_step);
    bool restart_from_checkpoint(int failed_step, int& target_step);
    bool reset_observables_after_restart(int failed_step);
    bool device_fields_finite(const char* label, int step);
    StepStatus run_step(bool screenprint, int step);
    void print_observables(const observableRecord& record);
    void release();
};
