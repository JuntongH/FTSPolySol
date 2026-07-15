#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cerrno>
#include <ctime>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cufft.h>

#include "errHandle.h"

using Complex = cufftDoubleComplex;

constexpr int ThreadsPerBlock = 256;
constexpr double Pi = 3.141592653589793238462643383279502884;
constexpr double UnsetSmearing = -1.0;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        errHandle::cuda(call);                                                 \
    } while (0)

#define CUFFT_CHECK(call)                                                      \
    do {                                                                       \
        errHandle::cufft(call);                                                \
    } while (0)

#define CUDA_LAUNCH(...)                                                       \
    do {                                                                       \
        __VA_ARGS__;                                                           \
        CUDA_CHECK(cudaGetLastError());                                        \
    } while (0)

inline int blocks(int n)
{
    return (n + ThreadsPerBlock - 1) / ThreadsPerBlock;
}

inline size_t linear_index(int x, int y, int z, int m2, int m3)
{
    return (static_cast<size_t>(x) * static_cast<size_t>(m2) + static_cast<size_t>(y)) *
           static_cast<size_t>(m3) + static_cast<size_t>(z);
}

inline size_t grid_size(int m1, int m2, int m3)
{
    return static_cast<size_t>(m1) * static_cast<size_t>(m2) * static_cast<size_t>(m3);
}

struct Segment {
    int type = 0;
    double charge = 0.0;
    double bond_length = 1.0;
};

struct IonSpecies {
    std::string name;
    double valence = 0.0;
    long long count = 0;
    double smearing_length = 0.0;
};

__host__ __device__ inline Complex make_complex(double x, double y)
{
    Complex z;
    z.x = x;
    z.y = y;
    return z;
}

__host__ __device__ inline Complex complex_zero()
{
    return make_complex(0.0, 0.0);
}

__host__ __device__ inline Complex complex_add(Complex a, Complex b)
{
    return make_complex(a.x + b.x, a.y + b.y);
}

__host__ __device__ inline Complex complex_sub(Complex a, Complex b)
{
    return make_complex(a.x - b.x, a.y - b.y);
}

__host__ __device__ inline Complex complex_mul(Complex a, Complex b)
{
    return make_complex(a.x * b.x - a.y * b.y,
                        a.x * b.y + a.y * b.x);
}

__host__ __device__ inline Complex complex_scale(Complex z, double scale)
{
    return make_complex(scale * z.x, scale * z.y);
}

__host__ __device__ inline double complex_abs2(Complex z)
{
    return z.x * z.x + z.y * z.y;
}

inline bool finite_complex(Complex z)
{
    return std::isfinite(z.x) && std::isfinite(z.y);
}

inline bool finite_real(double value)
{
    return std::isfinite(value);
}

__host__ __device__ inline Complex complex_inverse(Complex z)
{
    const double ax = fabs(z.x);
    const double ay = fabs(z.y);
    if (ax >= ay) {
        const double r = z.y / z.x;
        const double denom = z.x * (1.0 + r * r);
        return make_complex(1.0 / denom, -r / denom);
    }
    const double r = z.x / z.y;
    const double denom = z.y * (1.0 + r * r);
    return make_complex(r / denom, -1.0 / denom);
}

inline bool finite_complex_value(const char* name, Complex value, int step = -1)
{
    if (finite_complex(value)) return true;
    if (step >= 0) {
        std::fprintf(stderr, "error: non-finite value at step %d in %s = (% .10e, % .10e)\n",
                     step, name, value.x, value.y);
    } else {
        std::fprintf(stderr, "error: non-finite value in %s = (% .10e, % .10e)\n",
                     name, value.x, value.y);
    }
    return false;
}

inline bool finite_real_value(const char* name, double value, int step = -1)
{
    if (finite_real(value)) return true;
    if (step >= 0) {
        std::fprintf(stderr, "error: non-finite value at step %d in %s = % .10e\n", step, name, value);
    } else {
        std::fprintf(stderr, "error: non-finite value in %s = % .10e\n", name, value);
    }
    return false;
}

inline bool finite_complex_vector(const char* name, const std::vector<Complex>& values, int step = -1)
{
    for (size_t i = 0; i < values.size(); ++i) {
        const Complex value = values[i];
        if (!finite_complex(value)) {
            if (step >= 0) {
                std::fprintf(stderr, "error: non-finite value at step %d in %s[%zu] = (% .10e, % .10e)\n",
                             step, name, i, value.x, value.y);
            } else {
                std::fprintf(stderr, "error: non-finite value in %s[%zu] = (% .10e, % .10e)\n",
                             name, i, value.x, value.y);
            }
            return false;
        }
    }
    return true;
}

inline bool finite_complex_vectors(std::initializer_list<std::pair<const char*, const std::vector<Complex>*>> fields, int step = -1)
{
    for (const auto& field : fields) {
        if (!finite_complex_vector(field.first, *field.second, step)) return false;
    }
    return true;
}

__host__ __device__ inline Complex complex_exp(Complex z)
{
    const double r = exp(z.x);
    return make_complex(r * cos(z.y), r * sin(z.y));
}

inline Complex complex_neg_log(Complex z)
{
    const double r = std::sqrt(complex_abs2(z));
    return make_complex(-std::log(r), -std::atan2(z.y, z.x));
}

inline Complex complex_mean(Complex z, int size)
{
    return complex_scale(z, 1.0 / static_cast<double>(size));
}

inline Complex complex_lerp(Complex a, Complex b, double t)
{
    return make_complex(a.x + t * (b.x - a.x),
                        a.y + t * (b.y - a.y));
}

inline int periodic_index(int i, int n)
{
    return (i % n + n) % n;
}

inline int spectral_mode(int i, int m)
{
    return i > m / 2 ? i - m : i;
}

inline double wave_number(int i, int m, double length)
{
    return 2.0 * Pi * spectral_mode(i, m) / length;
}

inline bool line_content(char* line, char*& ptr)
{
    char* hash = std::strchr(line, '#');
    if (hash) *hash = '\0';
    ptr = line;
    while (std::isspace(static_cast<unsigned char>(*ptr))) ++ptr;
    return *ptr != '\0';
}

inline bool parse_mesh_header(const char* line, int& m1, int& m2, int& m3)
{
    char extra = 0;
    const int count = std::sscanf(line, "%d %d %d %c", &m1, &m2, &m3, &extra);
    return count == 3 && m1 > 0 && m2 > 0 && m3 > 0;
}

inline bool parse_box_lengths_header(const char* line, double& L1, double& L2, double& L3)
{
    char extra = 0;
    const int count = std::sscanf(line, "%lf %lf %lf %c", &L1, &L2, &L3, &extra);
    return count == 3 && L1 > 0.0 && L2 > 0.0 && L3 > 0.0;
}

inline bool write_mesh_header(FILE* fp, int m1, int m2, int m3, double L1, double L2, double L3)
{
    return std::fprintf(fp, "%d %d %d\n", m1, m2, m3) > 0 &&
           std::fprintf(fp, "%.17g %.17g %.17g\n", L1, L2, L3) > 0;
}

inline bool close_written_file(FILE* fp, const std::string& filename)
{
    if (std::fclose(fp) == 0) return true;
    return errHandle::message("cannot close " + filename);
}

inline bool replace_file(const std::string& source, const std::string& target)
{
    if (std::rename(source.c_str(), target.c_str()) == 0) return true;
    const std::string text = "cannot replace " + target + " with " + source + ": " + std::strerror(errno);
    std::remove(source.c_str());
    return errHandle::message(text);
}

inline bool valid_step_filename_template(const std::string& pattern)
{
    const size_t first = pattern.find("%d");
    if (first == std::string::npos || pattern.find("%d", first + 2) != std::string::npos) return false;
    const size_t slash = pattern.find_last_of("/\\");
    return slash == std::string::npos || first > slash;
}

inline std::string step_filename(const std::string& pattern, int step)
{
    const size_t pos = pattern.find("%d");
    if (pos == std::string::npos) return pattern;
    std::string out = pattern.substr(0, pos);
    out += std::to_string(step);
    out += pattern.substr(pos + 2);
    return out;
}

inline double gaussian_smearing(double length, double k2)
{
    return std::exp(-0.5 * length * length * k2);
}

inline double chain_propagator_kernel(double bond_length, double k2)
{
    return std::exp(-bond_length * bond_length * k2 / 6.0);
}

inline double langevin_noise_stddev(bool enabled, double lambda, double dt, double dV)
{
    return (enabled ? 1.0 : 0.0) * std::sqrt(2.0 * lambda * dt / dV);
}

inline double ideal_chemical_potential(long long molecule_count, double volume)
{
    return std::log(static_cast<double>(molecule_count) / volume);
}

inline double ion_self_chemical_potential(double bjerrum_length,
                                          double valence,
                                          double charge_smearing_length)
{
    return -bjerrum_length * valence * valence /
           (2.0 * std::sqrt(Pi) * charge_smearing_length);
}

inline void build_smearing_kernels(int size,
                                   const std::vector<IonSpecies>& species,
                                   const std::vector<double>& k2,
                                   std::vector<double>& kernels)
{
    kernels.assign(static_cast<size_t>(size) * species.size(), 0.0);
    for (size_t j = 0; j < species.size(); ++j) {
        double* out = kernels.data() + j * static_cast<size_t>(size);
        for (int i = 0; i < size; ++i) out[i] = gaussian_smearing(species[j].smearing_length, k2[i]);
    }
}

__device__ inline int thread_index()
{
    return blockIdx.x * blockDim.x + threadIdx.x;
}

__device__ inline int thread_stride()
{
    return blockDim.x * gridDim.x;
}

template <typename T>
inline void device_alloc(T*& ptr, size_t count)
{
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T)));
}

template <typename T>
inline void device_free(T*& ptr)
{
    if (ptr) {
        CUDA_CHECK(cudaFree(ptr));
        ptr = nullptr;
    }
}

template <typename T>
inline void copy_to_device(T* device, const T* host, size_t count)
{
    if (count == 0) return;
    CUDA_CHECK(cudaMemcpy(device, host, count * sizeof(T), cudaMemcpyHostToDevice));
}

template <typename T>
inline void copy_to_device(T* device, const std::vector<T>& host)
{
    copy_to_device(device, host.data(), host.size());
}

template <typename T>
inline void copy_to_host(std::vector<T>& host, const T* device)
{
    if (host.empty()) return;
    CUDA_CHECK(cudaMemcpy(host.data(), device, host.size() * sizeof(T), cudaMemcpyDeviceToHost));
}

inline bool write_one_complex_field_file(const std::string& filename, int m1, int m2, int m3, double L1, double L2, double L3, const Complex* field, size_t size)
{
    FILE* fp = std::fopen(filename.c_str(), "w");
    if (!fp) return errHandle::file(filename);
    if (!write_mesh_header(fp, m1, m2, m3, L1, L2, L3)) {
        std::fclose(fp);
        return errHandle::message("cannot write mesh header to " + filename);
    }
    for (size_t i = 0; i < size; ++i) {
        if (std::fprintf(fp, "%.10f  %.10f\n", field[i].x, field[i].y) < 0) {
            std::fclose(fp);
            return errHandle::message("cannot write " + filename);
        }
    }
    return close_written_file(fp, filename);
}

struct ComplexFieldTable {
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    bool header = false;
    std::vector<Complex> values;
};

inline bool read_complex_table_file(const std::string& filename,
                                    int default_m1,
                                    int default_m2,
                                    int default_m3,
                                    ComplexFieldTable& table)
{
    FILE* fp = std::fopen(filename.c_str(), "r");
    if (!fp) return errHandle::file(filename);

    table = ComplexFieldTable();
    table.m1 = default_m1;
    table.m2 = default_m2;
    table.m3 = default_m3;

    bool first_line = true;
    bool maybe_box_lengths = false;
    char line[4096];
    size_t line_number = 0;
    while (std::fgets(line, sizeof(line), fp)) {
        ++line_number;
        char* ptr = nullptr;
        if (!line_content(line, ptr)) continue;
        if (first_line) {
            int h1 = 0;
            int h2 = 0;
            int h3 = 0;
            if (parse_mesh_header(ptr, h1, h2, h3)) {
                table.m1 = h1;
                table.m2 = h2;
                table.m3 = h3;
                table.header = true;
                first_line = false;
                maybe_box_lengths = true;
                continue;
            }
            first_line = false;
        }
        if (maybe_box_lengths) {
            double L1 = 0.0;
            double L2 = 0.0;
            double L3 = 0.0;
            if (parse_box_lengths_header(ptr, L1, L2, L3)) {
                maybe_box_lengths = false;
                continue;
            }
            maybe_box_lengths = false;
        }

        double real_part = 0.0;
        double imag_part = 0.0;
        char extra = 0;
        const int count = std::sscanf(ptr, "%lf %lf %c", &real_part, &imag_part, &extra);
        if (count != 2) {
            std::fclose(fp);
            return errHandle::line(filename, line_number);
        }
        table.values.push_back(make_complex(real_part, imag_part));
    }
    std::fclose(fp);
    return true;
}

inline bool require_table_mesh(const ComplexFieldTable& table, const std::string& filename, int m1, int m2, int m3, const char* message)
{
    if (table.m1 == m1 && table.m2 == m2 && table.m3 == m3) return true;
    return errHandle::message(filename + ": " + message);
}

inline bool require_table_rows(const ComplexFieldTable& table, const std::string& filename, size_t expected)
{
    if (table.values.size() == expected) return true;
    return errHandle::rows(filename, table.values.size(), expected);
}

static __global__ void fillComplexKernel(Complex* data, Complex value, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) data[i] = value;
}

static __global__ void copyComplexKernel(Complex* dst, const Complex* src, int offset, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] = src[i + offset];
}

static __global__ void addComplexKernel(Complex* dst, const Complex* src, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] = complex_add(dst[i], src[i]);
}

static __global__ void addScaledComplexKernel(Complex* dst, const Complex* src, double scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] = complex_add(dst[i], complex_scale(src[i], scale));
}

static __global__ void addScaledRealKernel(Complex* dst, const double* src, double scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i].x += scale * src[i];
}

static __global__ void addComplexConstKernel(Complex* dst, Complex value, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] = complex_add(dst[i], value);
}

static __global__ void addRealConstKernel(double* dst, double value, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] += value;
}

static __global__ void addScaledComplexConstKernel(Complex* dst, const Complex* src, Complex scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) dst[i] = complex_add(dst[i], complex_mul(src[i], scale));
}

static __global__ void scaleRealKernel(Complex* data, const double* scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) data[i] = complex_scale(data[i], scale[i]);
}

static __global__ void scaleRealConstKernel(Complex* data, double scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) data[i] = complex_scale(data[i], scale);
}

static __global__ void scaleComplexKernel(Complex* data, const Complex* scale, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) data[i] = complex_mul(data[i], scale[i]);
}

static __global__ void expComplexKernel(Complex* out, const Complex* in, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) out[i] = complex_exp(in[i]);
}

static __global__ void expMinusComplexKernel(Complex* out, const Complex* in, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) out[i] = complex_exp(make_complex(-in[i].x, -in[i].y));
}

static __global__ void addITimesComplexKernel(Complex* dst, const Complex* src, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) {
        dst[i].x -= src[i].y;
        dst[i].y += src[i].x;
    }
}

static __global__ void minusIValenceBoltzmannKernel(Complex* out, const Complex* field, double valence, int size)
{
    for (int i = thread_index(); i < size; i += thread_stride()) {
        out[i] = complex_exp(make_complex(valence * field[i].y, -valence * field[i].x));
    }
}

static __global__ void normalizedDensityKernel(Complex* rho, const Complex* boltzmann, double density, Complex Q, int size)
{
    const Complex factor = complex_scale(complex_inverse(Q), density);
    for (int i = thread_index(); i < size; i += thread_stride()) rho[i] = complex_mul(boltzmann[i], factor);
}

static __global__ void accumulatePolymerDensityByChainDensityKernel(Complex* rho,
                                                                     const Complex* exp_W,
                                                                     const Complex* qF,
                                                                     const Complex* qB,
                                                                     double chain_density,
                                                                     Complex Q,
                                                                     double weight,
                                                                     int segment,
                                                                     int size)
{
    if (weight == 0.0) return;

    const Complex factor = complex_scale(complex_inverse(Q), chain_density * weight);

    for (int i = thread_index(); i < size; i += thread_stride()) {
        const int idx = i + segment * size;
        const Complex q_product = complex_mul(qF[idx], qB[idx]);
        const Complex weighted = complex_mul(q_product, exp_W[i]);
        rho[i] = complex_add(rho[i], complex_mul(weighted, factor));
    }
}

static __global__ void generateNoiseKernel(double* eta, double stddev, unsigned long long seed, int size)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= size) return;

    curandState state;
    curand_init(seed, i, 0, &state);
    eta[i] = stddev * curand_normal(&state);
}

static __global__ void finiteFieldArrayKernel(const Complex* fields, int field_count, int size, int* ok)
{
    for (int idx = thread_index(); idx < field_count * size; idx += thread_stride()) {
        const Complex z = fields[idx];
        if (!isfinite(z.x) || !isfinite(z.y)) atomicExch(ok, 0);
    }
}

template <unsigned int blockSize>
__device__ void warp_reduce_complex(volatile Complex* sdata, int tid)
{
    if (blockSize >= 64) {
        sdata[tid].x += sdata[tid + 32].x;
        sdata[tid].y += sdata[tid + 32].y;
    }
    if (blockSize >= 32) {
        sdata[tid].x += sdata[tid + 16].x;
        sdata[tid].y += sdata[tid + 16].y;
    }
    if (blockSize >= 16) {
        sdata[tid].x += sdata[tid + 8].x;
        sdata[tid].y += sdata[tid + 8].y;
    }
    if (blockSize >= 8) {
        sdata[tid].x += sdata[tid + 4].x;
        sdata[tid].y += sdata[tid + 4].y;
    }
    if (blockSize >= 4) {
        sdata[tid].x += sdata[tid + 2].x;
        sdata[tid].y += sdata[tid + 2].y;
    }
    if (blockSize >= 2) {
        sdata[tid].x += sdata[tid + 1].x;
        sdata[tid].y += sdata[tid + 1].y;
    }
}

template <unsigned int blockSize>
static __global__ void reduceComplexSumKernel(const Complex* in, Complex* out, unsigned int n)
{
    extern __shared__ Complex sdata[];

    const unsigned int tid = threadIdx.x;
    const unsigned int grid_size = blockSize * 2 * gridDim.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + tid;

    sdata[tid] = complex_zero();

    while (i < n) {
        Complex value = in[i];
        if (i + blockSize < n) value = complex_add(value, in[i + blockSize]);
        sdata[tid] = complex_add(sdata[tid], value);
        i += grid_size;
    }

    __syncthreads();

    if (blockSize >= 512) {
        if (tid < 256) sdata[tid] = complex_add(sdata[tid], sdata[tid + 256]);
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128) sdata[tid] = complex_add(sdata[tid], sdata[tid + 128]);
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64) sdata[tid] = complex_add(sdata[tid], sdata[tid + 64]);
        __syncthreads();
    }

    if (tid < 32) warp_reduce_complex<blockSize>(sdata, tid);
    if (tid == 0) out[blockIdx.x] = sdata[0];
}

template <unsigned int blockSize>
__device__ void warp_reduce_sum_double(volatile double* sdata, int tid)
{
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
    if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
    if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
    if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

template <unsigned int blockSize>
static __global__ void reduceDoubleSumKernel(const double* in, double* out, unsigned int n)
{
    extern __shared__ double sdata_sum[];

    const unsigned int tid = threadIdx.x;
    const unsigned int grid_size = blockSize * 2 * gridDim.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + tid;

    double local = 0.0;
    while (i < n) {
        local += in[i];
        if (i + blockSize < n) local += in[i + blockSize];
        i += grid_size;
    }

    sdata_sum[tid] = local;
    __syncthreads();

    if (blockSize >= 512) {
        if (tid < 256) sdata_sum[tid] += sdata_sum[tid + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128) sdata_sum[tid] += sdata_sum[tid + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64) sdata_sum[tid] += sdata_sum[tid + 64];
        __syncthreads();
    }

    if (tid < 32) warp_reduce_sum_double<blockSize>(sdata_sum, tid);
    if (tid == 0) out[blockIdx.x] = sdata_sum[0];
}

template <unsigned int blockSize>
__device__ void warp_reduce_max_double(volatile double* sdata, int tid)
{
    if (blockSize >= 64) sdata[tid] = fmax(sdata[tid], sdata[tid + 32]);
    if (blockSize >= 32) sdata[tid] = fmax(sdata[tid], sdata[tid + 16]);
    if (blockSize >= 16) sdata[tid] = fmax(sdata[tid], sdata[tid + 8]);
    if (blockSize >= 8) sdata[tid] = fmax(sdata[tid], sdata[tid + 4]);
    if (blockSize >= 4) sdata[tid] = fmax(sdata[tid], sdata[tid + 2]);
    if (blockSize >= 2) sdata[tid] = fmax(sdata[tid], sdata[tid + 1]);
}

template <unsigned int blockSize>
static __global__ void reduceComplexMaxAbs2Kernel(const Complex* in, double* out, unsigned int n)
{
    extern __shared__ double sdata_max[];

    const unsigned int tid = threadIdx.x;
    const unsigned int grid_size = blockSize * 2 * gridDim.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + tid;

    double local = 0.0;
    while (i < n) {
        const Complex a = in[i];
        double value = a.x * a.x + a.y * a.y;
        if (!isfinite(value)) value = INFINITY;
        local = fmax(local, value);

        if (i + blockSize < n) {
            const Complex b = in[i + blockSize];
            value = b.x * b.x + b.y * b.y;
            if (!isfinite(value)) value = INFINITY;
            local = fmax(local, value);
        }
        i += grid_size;
    }

    sdata_max[tid] = local;
    __syncthreads();

    if (blockSize >= 512) {
        if (tid < 256) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 256]);
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 128]);
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 64]);
        __syncthreads();
    }

    if (tid < 32) warp_reduce_max_double<blockSize>(sdata_max, tid);
    if (tid == 0) out[blockIdx.x] = sdata_max[0];
}


template <unsigned int blockSize>
static __global__ void reduceDoubleMaxKernel(const double* in, double* out, unsigned int n)
{
    extern __shared__ double sdata_max[];

    const unsigned int tid = threadIdx.x;
    const unsigned int grid_size = blockSize * 2 * gridDim.x;
    unsigned int i = blockIdx.x * (blockSize * 2) + tid;

    double local = 0.0;
    while (i < n) {
        double value = in[i];
        if (!isfinite(value)) value = INFINITY;
        local = fmax(local, value);

        if (i + blockSize < n) {
            value = in[i + blockSize];
            if (!isfinite(value)) value = INFINITY;
            local = fmax(local, value);
        }
        i += grid_size;
    }

    sdata_max[tid] = local;
    __syncthreads();

    if (blockSize >= 512) {
        if (tid < 256) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 256]);
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 128]);
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64) sdata_max[tid] = fmax(sdata_max[tid], sdata_max[tid + 64]);
        __syncthreads();
    }

    if (tid < 32) warp_reduce_max_double<blockSize>(sdata_max, tid);
    if (tid == 0) out[blockIdx.x] = sdata_max[0];
}

inline void fill_complex(Complex* data, Complex value, int size)
{
    CUDA_LAUNCH(fillComplexKernel<<<blocks(size), ThreadsPerBlock>>>(data, value, size));
}

inline void clear_complex(Complex* data, int size)
{
    fill_complex(data, complex_zero(), size);
}

inline void copy_complex_slice(Complex* dst, const Complex* src, int slice, int size)
{
    CUDA_LAUNCH(copyComplexKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, slice * size, size));
}

inline void copy_complex(Complex* dst, const Complex* src, int size)
{
    copy_complex_slice(dst, src, 0, size);
}

inline void add_complex(Complex* dst, const Complex* src, int size)
{
    CUDA_LAUNCH(addComplexKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, size));
}

inline void add_scaled_complex(Complex* dst, const Complex* src, double scale, int size)
{
    CUDA_LAUNCH(addScaledComplexKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, scale, size));
}

inline void add_scaled_i_complex(Complex* dst, const Complex* src, double scale, int size)
{
    CUDA_LAUNCH(addScaledComplexConstKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, make_complex(0.0, scale), size));
}

inline void add_scaled_complex_const(Complex* dst, const Complex* src, Complex scale, int size)
{
    CUDA_LAUNCH(addScaledComplexConstKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, scale, size));
}

inline void add_scaled_real(Complex* dst, const double* src, double scale, int size)
{
    CUDA_LAUNCH(addScaledRealKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, scale, size));
}

inline void add_complex_const(Complex* dst, Complex value, int size)
{
    CUDA_LAUNCH(addComplexConstKernel<<<blocks(size), ThreadsPerBlock>>>(dst, value, size));
}

inline void add_real_const(double* dst, double value, int size)
{
    CUDA_LAUNCH(addRealConstKernel<<<blocks(size), ThreadsPerBlock>>>(dst, value, size));
}

inline void add_i_times_complex(Complex* dst, const Complex* src, int size)
{
    CUDA_LAUNCH(addITimesComplexKernel<<<blocks(size), ThreadsPerBlock>>>(dst, src, size));
}

inline void scale_complex_real(Complex* data, const double* scale, int size)
{
    CUDA_LAUNCH(scaleRealKernel<<<blocks(size), ThreadsPerBlock>>>(data, scale, size));
}

inline void scale_complex_real_const(Complex* data, double scale, int size)
{
    CUDA_LAUNCH(scaleRealConstKernel<<<blocks(size), ThreadsPerBlock>>>(data, scale, size));
}

inline void multiply_complex(Complex* data, const Complex* scale, int size)
{
    CUDA_LAUNCH(scaleComplexKernel<<<blocks(size), ThreadsPerBlock>>>(data, scale, size));
}

inline void exp_complex(Complex* out, const Complex* in, int size)
{
    CUDA_LAUNCH(expComplexKernel<<<blocks(size), ThreadsPerBlock>>>(out, in, size));
}

inline void exp_minus_complex(Complex* out, const Complex* in, int size)
{
    CUDA_LAUNCH(expMinusComplexKernel<<<blocks(size), ThreadsPerBlock>>>(out, in, size));
}

inline void boltzmann_minus_i_valence(Complex* out, const Complex* field, double valence, int size)
{
    CUDA_LAUNCH(minusIValenceBoltzmannKernel<<<blocks(size), ThreadsPerBlock>>>(out, field, valence, size));
}

inline void normalized_density(Complex* rho, const Complex* boltzmann, double density, Complex Q, int size)
{
    CUDA_LAUNCH(normalizedDensityKernel<<<blocks(size), ThreadsPerBlock>>>(rho, boltzmann, density, Q, size));
}

inline void accumulate_polymer_density_by_chain_density(Complex* rho,
                                                        const Complex* exp_W,
                                                        const Complex* qF,
                                                        const Complex* qB,
                                                        double chain_density,
                                                        Complex Q,
                                                        double weight,
                                                        int segment,
                                                        int size)
{
    CUDA_LAUNCH(accumulatePolymerDensityByChainDensityKernel<<<blocks(size), ThreadsPerBlock>>>
                (rho, exp_W, qF, qB, chain_density, Q, weight, segment, size));
}

inline unsigned long long new_seed()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return static_cast<unsigned long long>(ts.tv_nsec ^ ts.tv_sec);
}

inline void generate_noise(double* eta, double stddev, int size)
{
    CUDA_LAUNCH(generateNoiseKernel<<<blocks(size), ThreadsPerBlock>>>(eta, stddev, new_seed(), size));
}

inline Complex sum_complex(const Complex* data, int size)
{
    if (size <= 0) return complex_zero();

    constexpr int threads = ThreadsPerBlock;
    const int needed_blocks = (size + threads * 2 - 1) / (threads * 2);
    const int block_count = std::min(needed_blocks, 256);

    Complex* partial = nullptr;
    Complex* result = nullptr;
    Complex host_result = complex_zero();

    device_alloc(partial, static_cast<size_t>(256));
    device_alloc(result, static_cast<size_t>(1));

    CUDA_LAUNCH(reduceComplexSumKernel<threads><<<block_count, threads, threads * sizeof(Complex)>>>
                (data, partial, static_cast<unsigned int>(size)));

    if (block_count > 1) {
        CUDA_LAUNCH(reduceComplexSumKernel<threads><<<1, threads, threads * sizeof(Complex)>>>
                    (partial, result, static_cast<unsigned int>(block_count)));
    } else {
        CUDA_CHECK(cudaMemcpy(result, partial, sizeof(Complex), cudaMemcpyDeviceToDevice));
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(&host_result, result, sizeof(Complex), cudaMemcpyDeviceToHost));
    device_free(partial);
    device_free(result);
    return host_result;
}

inline Complex mean_complex(const Complex* data, int size)
{
    return complex_mean(sum_complex(data, size), size);
}

inline double sum_real(const double* data, int size)
{
    if (size <= 0) return 0.0;

    constexpr int threads = ThreadsPerBlock;
    const int needed_blocks = (size + threads * 2 - 1) / (threads * 2);
    const int block_count = std::min(needed_blocks, 256);

    double* partial = nullptr;
    double* result = nullptr;
    double host_result = 0.0;

    device_alloc(partial, static_cast<size_t>(256));
    device_alloc(result, static_cast<size_t>(1));

    CUDA_LAUNCH(reduceDoubleSumKernel<threads><<<block_count, threads, threads * sizeof(double)>>>
                (data, partial, static_cast<unsigned int>(size)));

    if (block_count > 1) {
        CUDA_LAUNCH(reduceDoubleSumKernel<threads><<<1, threads, threads * sizeof(double)>>>
                    (partial, result, static_cast<unsigned int>(block_count)));
    } else {
        CUDA_CHECK(cudaMemcpy(result, partial, sizeof(double), cudaMemcpyDeviceToDevice));
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(&host_result, result, sizeof(double), cudaMemcpyDeviceToHost));
    device_free(partial);
    device_free(result);
    return host_result;
}

inline double mean_real(const double* data, int size)
{
    return size > 0 ? sum_real(data, size) / static_cast<double>(size) : 0.0;
}

inline void remove_complex_mean(Complex* data, int size)
{
    if (size <= 0) return;
    add_complex_const(data, complex_scale(mean_complex(data, size), -1.0), size);
}

inline void remove_real_mean(double* data, int size)
{
    if (size <= 0) return;
    add_real_const(data, -mean_real(data, size), size);
}

inline double max_complex_magnitude(const Complex* data, double scale, int size)
{
    if (size <= 0) return 0.0;

    constexpr int threads = ThreadsPerBlock;
    const int needed_blocks = (size + threads * 2 - 1) / (threads * 2);
    const int block_count = std::min(needed_blocks, 256);

    double* partial = nullptr;
    double* result = nullptr;
    double host_result = 0.0;

    device_alloc(partial, static_cast<size_t>(256));
    device_alloc(result, static_cast<size_t>(1));

    CUDA_LAUNCH(reduceComplexMaxAbs2Kernel<threads><<<block_count, threads, threads * sizeof(double)>>>
                (data, partial, static_cast<unsigned int>(size)));

    if (block_count > 1) {
        CUDA_LAUNCH(reduceDoubleMaxKernel<threads><<<1, threads, threads * sizeof(double)>>>
                    (partial, result, static_cast<unsigned int>(block_count)));
    } else {
        CUDA_CHECK(cudaMemcpy(result, partial, sizeof(double), cudaMemcpyDeviceToDevice));
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(&host_result, result, sizeof(double), cudaMemcpyDeviceToHost));
    device_free(partial);
    device_free(result);

    return std::fabs(scale) * std::sqrt(std::max(0.0, host_result));
}

inline void fft_plan_3d(cufftHandle* plan, int n1, int n2, int n3)
{
    CUFFT_CHECK(cufftPlan3d(plan, n1, n2, n3, CUFFT_Z2Z));
}

inline void fft_destroy(cufftHandle plan)
{
    CUFFT_CHECK(cufftDestroy(plan));
}

inline void fft_forward(cufftHandle plan, const Complex* src, Complex* dst)
{
    CUFFT_CHECK(cufftExecZ2Z(plan,
                             const_cast<Complex*>(src),
                             dst,
                             CUFFT_FORWARD));
}

inline void fft_inverse(cufftHandle plan, const Complex* src, Complex* dst)
{
    CUFFT_CHECK(cufftExecZ2Z(plan,
                             const_cast<Complex*>(src),
                             dst,
                             CUFFT_INVERSE));
}

inline void fft_apply_real_kernel(cufftHandle plan, const Complex* src, Complex* dst, const double* kernel, int size)
{
    fft_forward(plan, src, dst);
    scale_complex_real_const(dst, 1.0 / static_cast<double>(size), size);
    scale_complex_real(dst, kernel, size);
    fft_inverse(plan, dst, dst);
}

inline void fft_apply_real_kernel_scaled(cufftHandle plan, const Complex* src, Complex* dst, const double* kernel, double scale, int size)
{
    fft_apply_real_kernel(plan, src, dst, kernel, size);
    scale_complex_real_const(dst, scale, size);
}

inline void fft_smooth(cufftHandle plan, Complex* data, const double* kernel, int size)
{
    fft_apply_real_kernel(plan, data, data, kernel, size);
}

inline void copy_and_smooth(cufftHandle plan, Complex* dst, const Complex* src, const double* kernel, int size)
{
    copy_complex(dst, src, size);
    fft_smooth(plan, dst, kernel, size);
}

inline Complex ion_partition_function(cufftHandle plan,
                                      Complex* smoothed_psi,
                                      Complex* boltzmann,
                                      const Complex* psi,
                                      const double* kernel,
                                      double valence,
                                      int size)
{
    copy_and_smooth(plan, smoothed_psi, psi, kernel, size);
    boltzmann_minus_i_valence(boltzmann, smoothed_psi, valence, size);
    return mean_complex(boltzmann, size);
}

inline Complex ion_density_from_psi(cufftHandle plan,
                                    Complex* rho,
                                    Complex* smoothed_psi,
                                    Complex* boltzmann,
                                    const Complex* psi,
                                    const double* kernel,
                                    double valence,
                                    double density,
                                    int size)
{
    const Complex Q = ion_partition_function(plan, smoothed_psi, boltzmann, psi, kernel, valence, size);
    normalized_density(rho, boltzmann, density, Q, size);
    return Q;
}

inline bool finite_device_field_array(const Complex* fields, int field_count, int size, int* device_flag)
{
    if (field_count <= 0) return true;
    int host_flag = 1;
    CUDA_CHECK(cudaMemcpy(device_flag, &host_flag, sizeof(int), cudaMemcpyHostToDevice));
    CUDA_LAUNCH(finiteFieldArrayKernel<<<blocks(field_count * size), ThreadsPerBlock>>>(fields, field_count, size, device_flag));
    CUDA_CHECK(cudaMemcpy(&host_flag, device_flag, sizeof(int), cudaMemcpyDeviceToHost));
    return host_flag != 0;
}
