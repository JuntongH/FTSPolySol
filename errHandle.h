#pragma once

#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <cuda_runtime.h>
#include <cufft.h>

class errHandle {
public:
    static bool message(const std::string& text);
    static bool file(const std::string& path);
    static bool missing(const std::string& option);
    static bool unknown(const std::string& text);
    static bool invalid(const std::string& text);
    static bool parse(const std::string& text);
    static bool line(const std::string& file, size_t line);
    static bool rows(const std::string& file, size_t got, size_t expected);
    static bool mesh(const std::string& file, size_t rows);
    static void usage(const char* program);
    static void cuda(cudaError_t code);
    static void cufft(cufftResult code);
    static void launch(cudaError_t code);
};
