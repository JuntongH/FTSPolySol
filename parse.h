#pragma once

#include "common.h"

class parse {
public:
    struct options {
        const char* param_file = "param";
        std::string omega_file;
        std::string psi_file;
        std::string sequence_file;
        std::string run_id;
        std::string run_parent = ".";
        bool omega_override = false;
        bool psi_override = false;
        bool sequence_override = false;
        bool run_id_set = false;
        bool run_parent_set = false;
        bool screenprint = false;
    };

    static bool args(int argc, char** argv, options& out);
    static void usage(const char* program);
    static std::string trim(const std::string& text);
    static std::string strip_comment(const std::string& text);
    static std::string unquote(std::string value);
    static bool boolean(const std::string& value, bool& out);
    static bool real(const std::string& value, double& out);
    static bool integer64(const std::string& value, long long& out);
    static bool integer(const std::string& value, int& out);
    static bool string_value(const std::string& value, std::string& out);
    static bool counterion_number(const std::string& value, bool& automatic, long long& number);
};
