#include "parse.h"

bool parse::args(int argc, char** argv, options& out)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-p") == 0) {
            if (++i == argc) return errHandle::missing("-p");
            out.param_file = argv[i];
        } else if (std::strcmp(argv[i], "--omega") == 0) {
            if (++i == argc) return errHandle::missing("--omega");
            out.omega_file = argv[i];
            out.omega_override = true;
        } else if (std::strcmp(argv[i], "-psi") == 0 || std::strcmp(argv[i], "--psi") == 0) {
            if (++i == argc) return errHandle::missing("--psi");
            out.psi_file = argv[i];
            out.psi_override = true;
        } else if (std::strcmp(argv[i], "-i") == 0 || std::strcmp(argv[i], "--id") == 0) {
            if (++i == argc) return errHandle::missing("-i");
            out.run_id = argv[i];
            out.run_id_set = true;
        } else if (std::strcmp(argv[i], "-s") == 0 || std::strcmp(argv[i], "--dir") == 0 || std::strcmp(argv[i], "--run-dir") == 0) {
            if (++i == argc) return errHandle::missing("-s");
            out.run_parent = argv[i];
            out.run_parent_set = true;
        } else if (std::strcmp(argv[i], "--seq") == 0) {
            if (++i == argc) return errHandle::missing("--seq");
            out.sequence_file = argv[i];
            out.sequence_override = true;
        } else if (std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--charge") == 0 || std::strcmp(argv[i], "--charges") == 0 || std::strcmp(argv[i], "--sigma") == 0) {
            return errHandle::message("charge files are not used; use --seq");
        } else if (std::strcmp(argv[i], "--screenprint") == 0) {
            out.screenprint = true;
        } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            errHandle::unknown(std::string("option ") + argv[i]);
            usage(argv[0]);
            return false;
        }
    }
    return true;
}

void parse::usage(const char* program)
{
    errHandle::usage(program);
}

std::string parse::trim(const std::string& text)
{
    size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first]))) ++first;
    size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1]))) --last;
    return text.substr(first, last - first);
}

std::string parse::strip_comment(const std::string& text)
{
    bool quoted = false;
    char quote = '\0';
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if ((c == '"' || c == '\'') && (i == 0 || text[i - 1] != '\\')) {
            if (!quoted) {
                quoted = true;
                quote = c;
            } else if (quote == c) {
                quoted = false;
            }
        }
        if (c == '#' && !quoted) return text.substr(0, i);
    }
    return text;
}

std::string parse::unquote(std::string value)
{
    value = trim(value);
    if (value.size() >= 2) {
        const char a = value.front();
        const char b = value.back();
        if ((a == '"' && b == '"') || (a == '\'' && b == '\'')) value = value.substr(1, value.size() - 2);
    }
    return value;
}

bool parse::boolean(const std::string& value, bool& out)
{
    std::string text = unquote(value);
    for (char& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (text == "true" || text == "1" || text == "yes" || text == "on") {
        out = true;
        return true;
    }
    if (text == "false" || text == "0" || text == "no" || text == "off") {
        out = false;
        return true;
    }
    return false;
}

bool parse::real(const std::string& value, double& out)
{
    const std::string text = unquote(value);
    char* end = nullptr;
    out = std::strtod(text.c_str(), &end);
    return end != text.c_str() && *end == '\0';
}

bool parse::integer64(const std::string& value, long long& out)
{
    const std::string text = unquote(value);
    char* end = nullptr;
    out = std::strtoll(text.c_str(), &end, 10);
    return end != text.c_str() && *end == '\0';
}

bool parse::integer(const std::string& value, int& out)
{
    long long temp = 0;
    if (!integer64(value, temp)) return false;
    out = static_cast<int>(temp);
    return static_cast<long long>(out) == temp;
}

bool parse::string_value(const std::string& value, std::string& out)
{
    out = unquote(value);
    return true;
}

bool parse::counterion_number(const std::string& value, bool& automatic, long long& number)
{
    const std::string text = unquote(value);
    if (text == "auto_neutralize") {
        automatic = true;
        return true;
    }
    if (!integer64(text, number)) return false;
    automatic = false;
    return true;
}
