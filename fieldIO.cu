#include "fieldIO.h"

bool fieldIO::prepare(const parse::options& options, param& p)
{
    screenprint = options.screenprint;
    if (!options.run_id_set) {
        if (options.run_parent_set) {
            return errHandle::message("-s requires -i");
        }
        run.clear();
    } else {
        if (!valid_run_id(options.run_id)) {
            return errHandle::invalid("run id " + options.run_id);
        }
        run = join_path(options.run_parent, options.run_id);
        if (!make_directories(run)) return false;
        if (!copy_case_files(options, p, run)) return false;
    }
    input_omega = options.omega_override ? options.omega_file : p.omega_file;
    input_psi = options.psi_override ? options.psi_file : p.psi_file;
    zero_missing_default_omega = (!options.omega_override && input_omega == "omega.rf");
    source_mesh_known = false;
    source_m1 = source_m2 = source_m3 = 0;
    p.omega_file = run_output_path(run, p.omega_file);
    p.psi_file = run_output_path(run, p.psi_file);
    p.tmp_omega_file = run_output_path(run, p.tmp_omega_file);
    p.tmp_psi_file = run_output_path(run, p.tmp_psi_file);
    p.rho_type_file = run_output_path(run, p.rho_type_file);
    p.rhop_file = run_output_path(run, p.rhop_file);
    p.rhoc_file = run_output_path(run, p.rhoc_file);
    p.rhoc_poly_file = run_output_path(run, p.rhoc_poly_file);
    p.rhoc_ion_file = run_output_path(run, p.rhoc_ion_file);
    p.ion_density_prefix = run_output_path(run, p.ion_density_prefix);
    if (!p.step_omega_file.empty()) p.step_omega_file = run_output_path(run, p.step_omega_file);
    if (!p.step_psi_file.empty()) p.step_psi_file = run_output_path(run, p.step_psi_file);
    if (!p.step_rho_type_file.empty()) p.step_rho_type_file = run_output_path(run, p.step_rho_type_file);
    if (!p.step_rhop_file.empty()) p.step_rhop_file = run_output_path(run, p.step_rhop_file);
    if (!p.step_rhoc_file.empty()) p.step_rhoc_file = run_output_path(run, p.step_rhoc_file);
    if (!p.step_rhoc_poly_file.empty()) p.step_rhoc_poly_file = run_output_path(run, p.step_rhoc_poly_file);
    if (!p.step_rhoc_ion_file.empty()) p.step_rhoc_ion_file = run_output_path(run, p.step_rhoc_ion_file);
    if (!p.step_ion_density_prefix.empty()) p.step_ion_density_prefix = run_output_path(run, p.step_ion_density_prefix);
    if (!prepare_step_output_templates(p)) return false;
    p.observable_file = run_output_path(run, p.observable_file);
    p.drift_max_file = run_output_path(run, p.drift_max_file);
    return true;
}

const std::string& fieldIO::run_path() const
{
    return run;
}

const std::string& fieldIO::input_omega_file() const
{
    return input_omega;
}

const std::string& fieldIO::input_psi_file() const
{
    return input_psi;
}

bool fieldIO::path_exists_as_directory(const std::string& path)
{
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool fieldIO::make_directory(const std::string& path)
{
    if (path.empty()) return false;
    if (path_exists_as_directory(path)) return true;
    if (mkdir(path.c_str(), 0775) == 0) return true;
    if (errno == EEXIST && path_exists_as_directory(path)) return true;
    return errHandle::message("cannot create " + path);
}

bool fieldIO::make_directories(const std::string& path)
{
    if (path.empty()) return false;
    std::string current;
    size_t i = 0;
    if (path[0] == '/') {
        current = "/";
        i = 1;
    }
    while (i <= path.size()) {
        const size_t next = path.find('/', i);
        const std::string part = path.substr(i, next == std::string::npos ? std::string::npos : next - i);
        if (!part.empty() && part != ".") {
            if (!current.empty() && current.back() != '/') current += '/';
            current += part;
            if (!make_directory(current)) return false;
        }
        if (next == std::string::npos) break;
        i = next + 1;
    }
    return path_exists_as_directory(path) || make_directory(path);
}

bool fieldIO::copy_file(const std::string& source, const std::string& target)
{
    FILE* in = std::fopen(source.c_str(), "rb");
    if (!in) return errHandle::file(source);
    std::vector<char> data;
    char buffer[65536];
    while (true) {
        const size_t got = std::fread(buffer, 1, sizeof(buffer), in);
        data.insert(data.end(), buffer, buffer + got);
        if (got < sizeof(buffer)) {
            if (std::ferror(in)) {
                std::fclose(in);
                return errHandle::message("cannot read " + source);
            }
            break;
        }
    }
    std::fclose(in);
    FILE* out = std::fopen(target.c_str(), "wb");
    if (!out) return errHandle::file(target);
    if (!data.empty() && std::fwrite(data.data(), 1, data.size(), out) != data.size()) {
        std::fclose(out);
        return errHandle::message("cannot write " + target);
    }
    return close_written_file(out, target);
}

bool fieldIO::copy_case_files(const parse::options& options, const param& p, const std::string& run_path)
{
    if (!copy_file(options.param_file, join_path(run_path, options.run_id + ".param"))) return false;
    if (!p.sequence_file.empty() && !copy_file(p.sequence_file, join_path(run_path, options.run_id + ".sequence"))) return false;
    return true;
}

bool fieldIO::valid_run_id(const std::string& id)
{
    if (id.empty()) return false;
    if (id == "." || id == "..") return false;
    for (char c : id) if (c == '/' || c == '\\' || c == '\0') return false;
    return true;
}

std::string fieldIO::join_path(const std::string& a, const std::string& b)
{
    if (a.empty() || a == ".") return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

std::string fieldIO::directory_path(const std::string& path)
{
    const size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return std::string();
    if (pos == 0) return path.substr(0, 1);
    return path.substr(0, pos);
}

std::string fieldIO::basename_path(const std::string& path)
{
    const size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return path;
    return path.substr(pos + 1);
}

std::string fieldIO::sequence_directory_name(const std::string& pattern)
{
    const std::string base = basename_path(pattern);
    const size_t pos = base.find("%d");
    std::string name = pos == std::string::npos ? base : base.substr(0, pos);
    while (!name.empty() && (name.back() == '_' || name.back() == '-' || name.back() == '.')) {
        name.pop_back();
    }
    if (!name.empty()) return name;

    return "step_fields";
}

std::string fieldIO::run_output_path(const std::string& run_path, const std::string& name)
{
    if (run_path.empty()) return name;
    return join_path(run_path, basename_path(name));
}

bool fieldIO::prepare_step_output_template(std::string& pattern)
{
    if (pattern.empty() || !valid_step_filename_template(pattern)) return true;

    const std::string base = basename_path(pattern);
    const std::string parent = directory_path(pattern);
    const std::string folder = sequence_directory_name(base);
    const std::string full_folder = (!parent.empty() && basename_path(parent) == folder)
                                  ? parent
                                  : join_path(parent, folder);

    if (!make_directories(full_folder)) return false;
    pattern = join_path(full_folder, base);
    return true;
}

bool fieldIO::prepare_step_output_templates(param& p)
{
    if (p.Nso < 0) return true;
    return prepare_step_output_template(p.step_omega_file) &&
           prepare_step_output_template(p.step_psi_file) &&
           prepare_step_output_template(p.step_rho_type_file) &&
           prepare_step_output_template(p.step_rhop_file) &&
           prepare_step_output_template(p.step_rhoc_file) &&
           prepare_step_output_template(p.step_rhoc_poly_file) &&
           prepare_step_output_template(p.step_rhoc_ion_file) &&
           prepare_step_output_template(p.step_ion_density_prefix);
}

bool fieldIO::infer_mesh(size_t rows, const param& p, int& m1, int& m2, int& m3)
{
    const size_t target = grid_size(p.m1, p.m2, p.m3);
    if (rows == target) {
        m1 = p.m1;
        m2 = p.m2;
        m3 = p.m3;
        return true;
    }

    const double ratio = target > rows ? static_cast<double>(target) / rows : static_cast<double>(rows) / target;
    int scale = static_cast<int>(std::llround(std::cbrt(ratio)));
    if (scale > 1 && std::fabs(static_cast<double>(scale) * scale * scale - ratio) < 1.0e-8) {
        if (target > rows) {
            if (p.m1 % scale == 0 && p.m2 % scale == 0 && p.m3 % scale == 0) {
                const int a = p.m1 / scale;
                const int b = p.m2 / scale;
                const int c = p.m3 / scale;
                if (grid_size(a, b, c) == rows) {
                    m1 = a;
                    m2 = b;
                    m3 = c;
                    return true;
                }
            }
        } else {
            const int a = p.m1 * scale;
            const int b = p.m2 * scale;
            const int c = p.m3 * scale;
            if (grid_size(a, b, c) == rows) {
                m1 = a;
                m2 = b;
                m3 = c;
                return true;
            }
        }
    }

    const int cube = static_cast<int>(std::llround(std::cbrt(static_cast<double>(rows))));
    if (cube > 0 && grid_size(cube, cube, cube) == rows) {
        m1 = cube;
        m2 = cube;
        m3 = cube;
        return true;
    }
    return false;
}


namespace {

struct MultiComplexFieldTable {
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    bool header = false;
    int field_count = 0;
    std::vector<Complex> fields; // field-major: field * rows + row
};

bool parse_real_row(char* ptr, std::vector<double>& values)
{
    values.clear();
    char* p = ptr;
    while (*p != '\0') {
        while (std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (*p == '\0') break;
        char* end = nullptr;
        const double value = std::strtod(p, &end);
        if (end == p) return false;
        values.push_back(value);
        p = end;
    }
    return true;
}

bool read_multi_complex_table_file_any(const std::string& filename,
                                       int default_m1,
                                       int default_m2,
                                       int default_m3,
                                       MultiComplexFieldTable& table)
{
    FILE* fp = std::fopen(filename.c_str(), "r");
    if (!fp) return errHandle::file(filename);

    table = MultiComplexFieldTable();
    table.m1 = default_m1;
    table.m2 = default_m2;
    table.m3 = default_m3;

    bool first_line = true;
    bool maybe_box_lengths = false;
    std::vector<std::vector<Complex>> rows;
    char line[8192];
    size_t line_number = 0;
    while (std::fgets(line, sizeof(line), fp)) {
        ++line_number;
        char* ptr = nullptr;
        if (!line_content(line, ptr)) continue;
        if (first_line) {
            int h1 = 0, h2 = 0, h3 = 0;
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
            double L1 = 0.0, L2 = 0.0, L3 = 0.0;
            if (parse_box_lengths_header(ptr, L1, L2, L3)) {
                maybe_box_lengths = false;
                continue;
            }
            maybe_box_lengths = false;
        }

        std::vector<double> values;
        if (!parse_real_row(ptr, values)) { std::fclose(fp); return errHandle::line(filename, line_number); }
        if (values.empty() || values.size() % 2 != 0) {
            std::fclose(fp);
            return errHandle::message(filename + ": expected an even number of numeric columns per row");
        }
        const int row_field_count = static_cast<int>(values.size() / 2);
        if (table.field_count == 0) table.field_count = row_field_count;
        if (row_field_count != table.field_count) {
            std::fclose(fp);
            return errHandle::message(filename + ": inconsistent number of field columns");
        }
        std::vector<Complex> row(static_cast<size_t>(row_field_count));
        for (int f = 0; f < row_field_count; ++f) {
            row[static_cast<size_t>(f)] = make_complex(values[static_cast<size_t>(2*f)], values[static_cast<size_t>(2*f + 1)]);
        }
        rows.push_back(row);
    }
    std::fclose(fp);

    table.fields.assign(static_cast<size_t>(table.field_count) * rows.size(), complex_zero());
    for (size_t i = 0; i < rows.size(); ++i) {
        for (int f = 0; f < table.field_count; ++f) {
            table.fields[static_cast<size_t>(f) * rows.size() + i] = rows[i][static_cast<size_t>(f)];
        }
    }
    return true;
}

bool read_multi_complex_table_file(const std::string& filename,
                                   int default_m1,
                                   int default_m2,
                                   int default_m3,
                                   int expected_field_count,
                                   MultiComplexFieldTable& table)
{
    if (!read_multi_complex_table_file_any(filename, default_m1, default_m2, default_m3, table)) return false;
    if (table.field_count != expected_field_count) {
        return errHandle::message(filename + ": expected " + std::to_string(2 * expected_field_count) + " numeric columns per row");
    }
    return true;
}

size_t multi_complex_rows(const MultiComplexFieldTable& table)
{
    return table.field_count > 0 ? table.fields.size() / static_cast<size_t>(table.field_count) : 0;
}

bool resize_multi_complex_table(const MultiComplexFieldTable& in,
                                int m1,
                                int m2,
                                int m3,
                                MultiComplexFieldTable& out)
{
    if (in.field_count < 0) {
        return errHandle::message("invalid negative field count while resizing omega fields");
    }

    const size_t source_size = grid_size(in.m1, in.m2, in.m3);
    const size_t rows = multi_complex_rows(in);
    if (rows != source_size) {
        return errHandle::rows("omega field resize source", rows, source_size);
    }

    out = MultiComplexFieldTable();
    out.m1 = m1;
    out.m2 = m2;
    out.m3 = m3;
    out.header = true;
    out.field_count = in.field_count;

    const size_t target_size = grid_size(m1, m2, m3);
    out.fields.assign(static_cast<size_t>(std::max(0, in.field_count)) * target_size, complex_zero());

    for (int f = 0; f < in.field_count; ++f) {
        const size_t source_offset = static_cast<size_t>(f) * source_size;
        std::vector<Complex> source(in.fields.begin() + static_cast<ptrdiff_t>(source_offset),
                                    in.fields.begin() + static_cast<ptrdiff_t>(source_offset + source_size));
        std::vector<Complex> resized;
        trilinear::resize(source, in.m1, in.m2, in.m3, resized, m1, m2, m3);
        if (resized.size() != target_size) {
            return errHandle::message("omega field resize produced an unexpected number of grid points");
        }
        std::copy(resized.begin(), resized.end(),
                  out.fields.begin() + static_cast<ptrdiff_t>(static_cast<size_t>(f) * target_size));
    }

    return true;
}

int active_index_by_sign(const param& p, bool positive)
{
    for (int a = 0; a < p.omega_count(); ++a) {
        const double d = p.modes.active_eigenvalue(a);
        if (positive && d > 0.0) return a;
        if (!positive && d < 0.0) return a;
    }
    return -1;
}

bool map_legacy_ab_fields_to_active(const param& p, const MultiComplexFieldTable& table, std::vector<Complex>& omega)
{
    const int nfields = p.omega_count();
    const size_t size = static_cast<size_t>(p.size());
    omega.assign(static_cast<size_t>(std::max(0, nfields)) * size, complex_zero());
    if (nfields <= 0) return true;
    if (table.field_count != 1 && table.field_count != 2) {
        return errHandle::message("legacy A/B omega file must contain either one complex field or two complex fields");
    }
    if (multi_complex_rows(table) != size) {
        return errHandle::rows("legacy A/B omega field", multi_complex_rows(table), size);
    }

    const int composition_active = active_index_by_sign(p, false);
    const int density_active = active_index_by_sign(p, true);
    if (composition_active >= 0) {
        std::copy(table.fields.begin(), table.fields.begin() + static_cast<ptrdiff_t>(size),
                  omega.begin() + static_cast<ptrdiff_t>(composition_active) * static_cast<ptrdiff_t>(size));
    }
    if (table.field_count == 2 && density_active >= 0) {
        std::copy(table.fields.begin() + static_cast<ptrdiff_t>(size), table.fields.begin() + static_cast<ptrdiff_t>(2 * size),
                  omega.begin() + static_cast<ptrdiff_t>(density_active) * static_cast<ptrdiff_t>(size));
    }
    return true;
}

bool legacy_ab_fields_from_active(const param& p, const std::vector<Complex>& omega, std::vector<Complex>& legacy)
{
    const size_t size = static_cast<size_t>(p.size());
    legacy.assign(2 * size, complex_zero());
    const int composition_active = active_index_by_sign(p, false);
    const int density_active = active_index_by_sign(p, true);
    if (composition_active >= 0) {
        std::copy(omega.begin() + static_cast<ptrdiff_t>(composition_active) * static_cast<ptrdiff_t>(size),
                  omega.begin() + static_cast<ptrdiff_t>(composition_active + 1) * static_cast<ptrdiff_t>(size),
                  legacy.begin());
    }
    if (density_active >= 0) {
        std::copy(omega.begin() + static_cast<ptrdiff_t>(density_active) * static_cast<ptrdiff_t>(size),
                  omega.begin() + static_cast<ptrdiff_t>(density_active + 1) * static_cast<ptrdiff_t>(size),
                  legacy.begin() + static_cast<ptrdiff_t>(size));
    }
    return true;
}


bool write_multi_complex_fields_file(const std::string& filename,
                                     int m1,
                                     int m2,
                                     int m3,
                                     double L1,
                                     double L2,
                                     double L3,
                                     const Complex* fields,
                                     int field_count,
                                     size_t size)
{
    FILE* fp = std::fopen(filename.c_str(), "w");
    if (!fp) return errHandle::file(filename);
    if (!write_mesh_header(fp, m1, m2, m3, L1, L2, L3)) {
        std::fclose(fp);
        return errHandle::message("cannot write mesh header to " + filename);
    }
    for (size_t i = 0; i < size; ++i) {
        for (int f = 0; f < field_count; ++f) {
            const Complex z = fields[static_cast<size_t>(f) * size + i];
            if (std::fprintf(fp, "%s%.10f  %.10f", f == 0 ? "" : "  ", z.x, z.y) < 0) {
                std::fclose(fp);
                return errHandle::message("cannot write " + filename);
            }
        }
        if (std::fprintf(fp, "\n") < 0) {
            std::fclose(fp);
            return errHandle::message("cannot write " + filename);
        }
    }
    return close_written_file(fp, filename);
}

} // namespace

bool fieldIO::read_omega_fields(const param& p, std::vector<Complex>& omega) const
{
    const int nfields = p.omega_count();
    const int size = p.size();
    omega.assign(static_cast<size_t>(std::max(0, nfields)) * static_cast<size_t>(size), complex_zero());
    if (nfields <= 0) return true;

    errno = 0;
    FILE* probe = std::fopen(input_omega.c_str(), "r");
    if (!probe) {
        if (zero_missing_default_omega && errno == ENOENT) {
            if (screenprint) {
                std::printf("%s not found; initializing omega fields to zero.\n", input_omega.c_str());
            }
            return true;
        }
        return errHandle::file(input_omega);
    }
    std::fclose(probe);

    MultiComplexFieldTable table;
    if (p.legacy_ab_field_order()) {
        if (!read_multi_complex_table_file_any(input_omega, p.m1, p.m2, p.m3, table)) return false;
    } else {
        if (!read_multi_complex_table_file(input_omega, p.m1, p.m2, p.m3, nfields, table)) return false;
    }

    size_t expected = grid_size(table.m1, table.m2, table.m3);
    const size_t rows = multi_complex_rows(table);
    if (!table.header && rows != expected) {
        if (!infer_mesh(rows, p, table.m1, table.m2, table.m3)) return errHandle::mesh(input_omega, rows);
        expected = grid_size(table.m1, table.m2, table.m3);
        if (screenprint) std::printf("inferred mesh of %s as %d %d %d\n", input_omega.c_str(), table.m1, table.m2, table.m3);
    }
    if (rows != expected) return errHandle::rows(input_omega, rows, expected);

    source_mesh_known = true;
    source_m1 = table.m1;
    source_m2 = table.m2;
    source_m3 = table.m3;

    MultiComplexFieldTable final_table = table;
    if (table.m1 != p.m1 || table.m2 != p.m2 || table.m3 != p.m3) {
        if (screenprint) std::printf("resizing %s from %d %d %d to %d %d %d\n", input_omega.c_str(), table.m1, table.m2, table.m3, p.m1, p.m2, p.m3);
        if (!resize_multi_complex_table(table, p.m1, p.m2, p.m3, final_table)) return false;
    }

    if (p.legacy_ab_field_order()) {
        if (final_table.field_count == 1 && screenprint) {
            std::printf("%s has one complex field; reading it as the legacy A/B composition field and setting the legacy density field to zero.\n", input_omega.c_str());
        }
        return map_legacy_ab_fields_to_active(p, final_table, omega);
    }

    omega.swap(final_table.fields);
    return true;
}

bool fieldIO::read_psi_field(const param& p, std::vector<Complex>& psi) const
{
    FILE* probe = std::fopen(input_psi.c_str(), "r");
    if (!probe) {
        if (screenprint) std::printf("%s not found; initializing Psi to zero.\n", input_psi.c_str());
        for (Complex& z : psi) z = complex_zero();
        return true;
    }
    std::fclose(probe);

    ComplexFieldTable table;
    if (!read_complex_table_file(input_psi, p.m1, p.m2, p.m3, table)) return false;

    size_t expected = grid_size(table.m1, table.m2, table.m3);
    if (table.values.size() != expected && !table.header && source_mesh_known) {
        const size_t source_expected = grid_size(source_m1, source_m2, source_m3);
        if (table.values.size() == source_expected) {
            table.m1 = source_m1;
            table.m2 = source_m2;
            table.m3 = source_m3;
            expected = source_expected;
        }
    }
    if (table.values.size() != expected && !table.header) {
        if (!infer_mesh(table.values.size(), p, table.m1, table.m2, table.m3)) return errHandle::mesh(input_psi, table.values.size());
        expected = grid_size(table.m1, table.m2, table.m3);
        if (screenprint) std::printf("inferred mesh of %s as %d %d %d\n", input_psi.c_str(), table.m1, table.m2, table.m3);
    }
    if (table.values.size() != expected) return errHandle::rows(input_psi, table.values.size(), expected);

    if (table.m1 != p.m1 || table.m2 != p.m2 || table.m3 != p.m3) {
        if (screenprint) std::printf("resizing %s from %d %d %d to %d %d %d\n", input_psi.c_str(), table.m1, table.m2, table.m3, p.m1, p.m2, p.m3);
        trilinear::resize(table.values, table.m1, table.m2, table.m3, psi, p.m1, p.m2, p.m3);
    } else {
        psi = table.values;
    }
    return true;
}

bool fieldIO::read_omega_fields_exact(const param& p, const std::string& filename, std::vector<Complex>& omega)
{
    const int nfields = p.omega_count();
    const size_t expected = grid_size(p.m1, p.m2, p.m3);
    omega.assign(static_cast<size_t>(std::max(0, nfields)) * expected, complex_zero());
    if (nfields <= 0) return true;

    MultiComplexFieldTable table;
    if (p.legacy_ab_field_order()) {
        if (!read_multi_complex_table_file_any(filename, p.m1, p.m2, p.m3, table)) return false;
    } else {
        if (!read_multi_complex_table_file(filename, p.m1, p.m2, p.m3, nfields, table)) return false;
    }
    const size_t rows = multi_complex_rows(table);
    if (table.m1 != p.m1 || table.m2 != p.m2 || table.m3 != p.m3) return errHandle::message(filename + ": checkpoint mesh does not match current mesh");
    if (rows != expected) return errHandle::rows(filename, rows, expected);
    if (p.legacy_ab_field_order()) return map_legacy_ab_fields_to_active(p, table, omega);
    omega.swap(table.fields);
    return true;
}

bool fieldIO::read_complex_field_exact(const param& p, const std::string& filename, std::vector<Complex>& field)
{
    ComplexFieldTable table;
    if (!read_complex_table_file(filename, p.m1, p.m2, p.m3, table)) return false;
    const size_t expected = grid_size(p.m1, p.m2, p.m3);
    if (!require_table_mesh(table, filename, p.m1, p.m2, p.m3, "checkpoint mesh does not match current mesh")) return false;
    if (!require_table_rows(table, filename, expected)) return false;
    field.swap(table.values);
    return true;
}

bool fieldIO::write_omega_fields(const param& p, const std::string& filename, const std::vector<Complex>& omega)
{
    if (p.omega_count() <= 0) return write_multi_complex_fields_file(filename, p.m1, p.m2, p.m3, p.L1, p.L2, p.L3, nullptr, 0, static_cast<size_t>(p.size()));
    if (p.legacy_ab_field_order()) {
        std::vector<Complex> legacy;
        if (!legacy_ab_fields_from_active(p, omega, legacy)) return false;
        return write_multi_complex_fields_file(filename, p.m1, p.m2, p.m3, p.L1, p.L2, p.L3, legacy.data(), 2, static_cast<size_t>(p.size()));
    }
    return write_multi_complex_fields_file(filename, p.m1, p.m2, p.m3, p.L1, p.L2, p.L3, omega.data(), p.omega_count(), static_cast<size_t>(p.size()));
}

bool fieldIO::write_type_fields(const param& p, const std::string& filename, const std::vector<Complex>& fields, int field_count)
{
    return write_multi_complex_fields_file(filename, p.m1, p.m2, p.m3, p.L1, p.L2, p.L3, fields.data(), field_count, static_cast<size_t>(p.size()));
}

bool fieldIO::write_complex_field(const param& p, const std::string& filename, const Complex* field, size_t size)
{
    return write_one_complex_field_file(filename, p.m1, p.m2, p.m3, p.L1, p.L2, p.L3, field, size);
}

bool fieldIO::write_complex_field(const param& p, const std::string& filename, const std::vector<Complex>& field)
{
    return write_complex_field(p, filename, field.data(), field.size());
}

bool fieldIO::write_ion_densities_with_prefix(const param& p, const ion& ions, const std::string& prefix, const std::vector<Complex>& ion_rho)
{
    bool ok = true;
    const int size = p.size();
    const std::vector<IonSpecies>& list = ions.all();
    for (size_t j = 0; j < list.size(); ++j) {
        char filename[1024];
        std::snprintf(filename, sizeof(filename), "%s_%zu_%s.rf", prefix.c_str(), j, list[j].name.c_str());
        const Complex* begin = ion_rho.data() + j * static_cast<size_t>(size);
        ok = write_complex_field(p, std::string(filename), begin, size) && ok;
    }
    return ok;
}

bool fieldIO::write_ion_densities(const param& p, const ion& ions, const std::vector<Complex>& ion_rho)
{
    return write_ion_densities_with_prefix(p, ions, p.ion_density_prefix, ion_rho);
}

bool fieldIO::step_outputs_enabled(const param& p)
{
    return p.Nso >= 0 &&
           (!p.step_omega_file.empty() ||
            !p.step_psi_file.empty() ||
            !p.step_rho_type_file.empty() ||
            !p.step_rhop_file.empty() ||
            !p.step_rhoc_file.empty() ||
            !p.step_rhoc_poly_file.empty() ||
            !p.step_rhoc_ion_file.empty() ||
            !p.step_ion_density_prefix.empty());
}

bool fieldIO::write_step_outputs(const param& p, const ion& ions, int output_index, const std::vector<Complex>& omega, const std::vector<Complex>& psi, const std::vector<Complex>& rho_type, const std::vector<Complex>& rho_p, const std::vector<Complex>& rho_c, const std::vector<Complex>& rho_c_poly, const std::vector<Complex>& rho_c_ion, const std::vector<Complex>& ion_rho)
{
    bool ok = true;
    if (!p.step_omega_file.empty()) ok = write_omega_fields(p, step_filename(p.step_omega_file, output_index), omega) && ok;
    if (!p.step_psi_file.empty()) ok = write_complex_field(p, step_filename(p.step_psi_file, output_index), psi) && ok;
    if (!p.step_rho_type_file.empty()) ok = write_type_fields(p, step_filename(p.step_rho_type_file, output_index), rho_type, p.monomer_count()) && ok;
    if (!p.step_rhop_file.empty()) ok = write_complex_field(p, step_filename(p.step_rhop_file, output_index), rho_p) && ok;
    if (!p.step_rhoc_file.empty()) ok = write_complex_field(p, step_filename(p.step_rhoc_file, output_index), rho_c) && ok;
    if (!p.step_rhoc_poly_file.empty()) ok = write_complex_field(p, step_filename(p.step_rhoc_poly_file, output_index), rho_c_poly) && ok;
    if (!p.step_rhoc_ion_file.empty()) ok = write_complex_field(p, step_filename(p.step_rhoc_ion_file, output_index), rho_c_ion) && ok;
    if (!p.step_ion_density_prefix.empty()) ok = write_ion_densities_with_prefix(p, ions, step_filename(p.step_ion_density_prefix, output_index), ion_rho) && ok;
    return ok;
}

bool fieldIO::write_all(const param& p, const ion& ions, const std::vector<Complex>& omega, const std::vector<Complex>& psi, const std::vector<Complex>& rho_type, const std::vector<Complex>& rho_p, const std::vector<Complex>& rho_c, const std::vector<Complex>& rho_c_poly, const std::vector<Complex>& rho_c_ion, const std::vector<Complex>& ion_rho) const
{
    return write_omega_fields(p, p.omega_file, omega) &&
           write_complex_field(p, p.psi_file, psi) &&
           write_type_fields(p, p.rho_type_file, rho_type, p.monomer_count()) &&
           write_complex_field(p, p.rhop_file, rho_p) &&
           write_complex_field(p, p.rhoc_file, rho_c) &&
           write_complex_field(p, p.rhoc_poly_file, rho_c_poly) &&
           write_complex_field(p, p.rhoc_ion_file, rho_c_ion) &&
           write_ion_densities(p, ions, ion_rho);
}
