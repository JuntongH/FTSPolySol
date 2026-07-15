#pragma once

#include "param.h"
#include "deviceState.h"

class ion {
public:
    bool build(const param& p, double total_polymer_charge);
    void build_kernels(const param& p, const std::vector<double>& k2, std::vector<double>& ion_kernels) const;
    void compute_densities(cufftHandle plan, deviceState& d, const param& p, const Complex* psi) const;
    void print(const param& p) const;
    const std::vector<IonSpecies>& all() const;
    bool empty() const;
    size_t count() const;

private:
    std::vector<IonSpecies> species;

    static bool nearly_equal(double a, double b, double tol);
    static bool integer_from_double(double value, long long& out);
    static std::string clean_name(std::string name);
    void append(const std::string& name, double valence, long long count, double smearing_length, bool combine);
};
