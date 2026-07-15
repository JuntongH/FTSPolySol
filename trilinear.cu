#include "trilinear.h"

Complex trilinear::value(const std::vector<Complex>& in, int x, int y, int z, int n1, int n2, int n3)
{
    x = periodic_index(x, n1);
    y = periodic_index(y, n2);
    z = periodic_index(z, n3);
    return in[static_cast<size_t>(linear_index(x, y, z, n2, n3))];
}

void trilinear::resize(const std::vector<Complex>& in,
                       int n1,
                       int n2,
                       int n3,
                       std::vector<Complex>& out,
                       int m1,
                       int m2,
                       int m3)
{
    if (n1 == m1 && n2 == m2 && n3 == m3) {
        out = in;
        return;
    }

    out.assign(grid_size(m1, m2, m3), complex_zero());

    for (int x = 0; x < m1; ++x) {
        const double sx = static_cast<double>(x) * n1 / m1;
        const int x0 = static_cast<int>(std::floor(sx));
        const int x1 = n1 == 1 ? 0 : x0 + 1;
        const double tx = n1 == 1 ? 0.0 : sx - x0;

        for (int y = 0; y < m2; ++y) {
            const double sy = static_cast<double>(y) * n2 / m2;
            const int y0 = static_cast<int>(std::floor(sy));
            const int y1 = n2 == 1 ? 0 : y0 + 1;
            const double ty = n2 == 1 ? 0.0 : sy - y0;

            for (int z = 0; z < m3; ++z) {
                const double sz = static_cast<double>(z) * n3 / m3;
                const int z0 = static_cast<int>(std::floor(sz));
                const int z1 = n3 == 1 ? 0 : z0 + 1;
                const double tz = n3 == 1 ? 0.0 : sz - z0;

                const Complex c000 = value(in, x0, y0, z0, n1, n2, n3);
                const Complex c001 = value(in, x0, y0, z1, n1, n2, n3);
                const Complex c010 = value(in, x0, y1, z0, n1, n2, n3);
                const Complex c011 = value(in, x0, y1, z1, n1, n2, n3);
                const Complex c100 = value(in, x1, y0, z0, n1, n2, n3);
                const Complex c101 = value(in, x1, y0, z1, n1, n2, n3);
                const Complex c110 = value(in, x1, y1, z0, n1, n2, n3);
                const Complex c111 = value(in, x1, y1, z1, n1, n2, n3);

                const Complex c00 = complex_lerp(c000, c001, tz);
                const Complex c01 = complex_lerp(c010, c011, tz);
                const Complex c10 = complex_lerp(c100, c101, tz);
                const Complex c11 = complex_lerp(c110, c111, tz);
                const Complex c0 = complex_lerp(c00, c01, ty);
                const Complex c1 = complex_lerp(c10, c11, ty);
                out[static_cast<size_t>(linear_index(x, y, z, m2, m3))] = complex_lerp(c0, c1, tx);
            }
        }
    }
}
