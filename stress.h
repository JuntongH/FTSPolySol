#pragma once

#include "polymer.h"
#include "ion.h"

class stress {
public:
    enum Component {
        XX = 0,
        YY = 1,
        ZZ = 2,
        XY = 3,
        XZ = 4,
        YZ = 5,
        NumComponents = 6
    };

    stress(const param& p, const polymerSpecies& poly, const ion& ions);
    ~stress();

    bool initialize();
    void release();

    bool enabled() const;
    size_t count() const;
    std::vector<std::string> labels() const;

    std::vector<Complex> measure(cufftHandle plan,
                                 deviceState& d,
                                 const std::vector<Complex>& Qp,
                                 const Complex* omega,
                                 const Complex* psi) const;

private:
    const param& p;
    const polymerSpecies& poly;
    const ion& ions;

    std::vector<Component> active_components;
    double mode_term[NumComponents] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double* kab = nullptr;
    double* minus_kab = nullptr;

    mutable Complex* scratch_a = nullptr;
    mutable Complex* scratch_b = nullptr;
    mutable Complex* scratch_c = nullptr;
    mutable Complex* dfield = nullptr;
    mutable Complex* dpsi = nullptr;

    void build_active_components();
    bool build_kernels();
    bool allocate_scratch();

    const double* kab_component(Component c) const;
    const double* minus_kab_component(Component c) const;


    Complex ideal_and_measure_term(Component c) const;
    Complex electrostatic_field_term(cufftHandle plan,
                                     deviceState& d,
                                     const Complex* psi,
                                     Component c) const;
    Complex chain_connectivity_term(cufftHandle plan,
                                    deviceState& d,
                                    const polymerSpecies::Species& sp,
                                    Complex Qp,
                                    Component c) const;
    Complex polymer_field_coupling_term(cufftHandle plan,
                                        deviceState& d,
                                        const polymerSpecies::Species& sp,
                                        Complex Qp,
                                        Component c) const;
    Complex ion_field_coupling_term(cufftHandle plan,
                                    deviceState& d,
                                    const Complex* psi,
                                    Component c) const;
};
