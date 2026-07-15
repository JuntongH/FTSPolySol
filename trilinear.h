#pragma once

#include "common.h"

class trilinear {
public:
    static void resize(const std::vector<Complex>& in,
                       int n1,
                       int n2,
                       int n3,
                       std::vector<Complex>& out,
                       int m1,
                       int m2,
                       int m3);

private:
    static Complex value(const std::vector<Complex>& in, int x, int y, int z, int n1, int n2, int n3);
};
