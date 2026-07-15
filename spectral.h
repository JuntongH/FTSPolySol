#pragma once

#include "polymer.h"

class spectral {
public:
    static void build_kernels(const param& p, const polymerSpecies& poly, std::vector<double>& Gamma, std::vector<double>& PHI, std::vector<double>& k2);
};
