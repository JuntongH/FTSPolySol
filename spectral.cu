#include "spectral.h"

void spectral::build_kernels(const param& p,
                             const polymerSpecies& poly,
                             std::vector<double>& Gamma,
                             std::vector<double>& PHI,
                             std::vector<double>& k2)
{
    const int size = p.size();
    const size_t bond_count = std::max<size_t>(1, poly.bond_kernel_count());
    Gamma.resize(static_cast<size_t>(size));
    PHI.assign(bond_count * static_cast<size_t>(size), 1.0);
    k2.resize(static_cast<size_t>(size));

    for (int x = 0; x < p.m1; ++x) {
        const double kx = wave_number(x, p.m1, p.L1);
        for (int y = 0; y < p.m2; ++y) {
            const double ky = wave_number(y, p.m2, p.L2);
            for (int z = 0; z < p.m3; ++z) {
                const double kz = wave_number(z, p.m3, p.L3);
                const int i = static_cast<int>(linear_index(x, y, z, p.m2, p.m3));
                const double ksq = kx * kx + ky * ky + kz * kz;
                k2[static_cast<size_t>(i)] = ksq;
                Gamma[static_cast<size_t>(i)] = gaussian_smearing(p.a, ksq);

                for (const polymerSpecies::Species& sp : poly.all()) {
                    for (int s = 0; s < sp.N - 1; ++s) {
                        const size_t offset = (sp.bond_offset + static_cast<size_t>(s)) * static_cast<size_t>(size) + static_cast<size_t>(i);
                        PHI[offset] = chain_propagator_kernel(sp.bond_length(s), ksq);
                    }
                }
            }
        }
    }
}
