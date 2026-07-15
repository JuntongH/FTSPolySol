#include "errHandle.h"

namespace {

bool emit(const std::string& text)
{
    std::fprintf(stderr, "error: %s\n", text.c_str());
    return false;
}

}

bool errHandle::message(const std::string& text)
{
    return emit(text);
}

bool errHandle::file(const std::string& path)
{
    std::string text = "cannot open " + path;
    if (errno) text += ": " + std::string(std::strerror(errno));
    return emit(text);
}

bool errHandle::missing(const std::string& option)
{
    return emit("missing value for " + option);
}

bool errHandle::unknown(const std::string& text)
{
    return emit("unknown " + text);
}

bool errHandle::invalid(const std::string& text)
{
    return emit("invalid " + text);
}

bool errHandle::parse(const std::string& text)
{
    return emit("cannot parse " + text);
}

bool errHandle::line(const std::string& file, size_t line_number)
{
    return emit("bad line " + std::to_string(line_number) + " in " + file);
}

bool errHandle::rows(const std::string& file, size_t got, size_t expected)
{
    return emit(file + ": rows " + std::to_string(got) + ", expected " + std::to_string(expected));
}

bool errHandle::mesh(const std::string& file, size_t rows)
{
    return emit(file + ": cannot infer mesh from " + std::to_string(rows) + " rows");
}

void errHandle::usage(const char* program)
{
    std::fprintf(stderr, "usage: %s [-p param] [-i id] [-s dir] [--omega omega_file] [--psi psi_file] [--seq sequence_file] [--screenprint]\n", program);
}

void errHandle::cuda(cudaError_t code)
{
    if (code == cudaSuccess) return;
    emit("cuda: " + std::string(cudaGetErrorString(code)));
    std::exit(EXIT_FAILURE);
}

void errHandle::cufft(cufftResult code)
{
    if (code == CUFFT_SUCCESS) return;
    emit("cufft: code " + std::to_string(static_cast<int>(code)));
    std::exit(EXIT_FAILURE);
}

void errHandle::launch(cudaError_t code)
{
    if (code == cudaSuccess) return;
    emit("kernel: " + std::string(cudaGetErrorString(code)));
    std::exit(EXIT_FAILURE);
}
