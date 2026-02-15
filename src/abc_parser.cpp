#include <set>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cmath>

#include "logging.h"
#include "file.h"
#include "abc_parser.h"

// ============================================================================
// Helper utilities
// ============================================================================

static std::string str_trim(const std::string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static std::string strip_comment(const std::string& line)
{
    size_t pos = line.find("//");
    if (pos == std::string::npos)
        return line;
    std::string result = line.substr(0, pos);
    while (!result.empty() && std::isspace((unsigned char)result.back()))
        result.pop_back();
    return result;
}

static bool is_skip_line(const std::string& line)
{
    if (line.empty())
        return true;
    size_t first_non_ws = line.find_first_not_of(" \t");
    if (first_non_ws == std::string::npos)
        return true;
    if (line.size() > first_non_ws + 1 && line[first_non_ws] == '/' && line[first_non_ws + 1] == '/')
        return true;
    return false;
}

static bool has_prefix_ci(const std::string& s, const char* prefix)
{
    size_t len = std::strlen(prefix);
    if (s.size() < len)
        return false;
    for (size_t i = 0; i < len; i++)
    {
        if (std::tolower((unsigned char)s[i]) != std::tolower((unsigned char)prefix[i]))
            return false;
    }
    return true;
}

FileFormat detect_format(const std::string& first_line)
{
    if (has_prefix_ci(first_line, "ABCD "))
        return FORMAT_ABCD;
    if (has_prefix_ci(first_line, "ABC2 "))
        return FORMAT_ABC2;
    if (has_prefix_ci(first_line, "ABC "))
        return FORMAT_ABC;
    return FORMAT_UNKNOWN;
}

// ============================================================================
// Prime sieve for ABC2 "primes from X to Y" support
// ============================================================================

static std::vector<int64_t> sieve_primes(int64_t min_val, int64_t max_val)
{
    std::vector<int64_t> result;
    if (max_val < 2) return result;
    if (min_val < 2) min_val = 2;
    int64_t range = max_val - min_val + 1;
    if (range <= 0) return result;
    if (range > MAX_SIEVE_RANGE)
        return result; // Range too large — caller must check and warn

    int64_t sqrt_max = (int64_t)std::sqrt((double)max_val) + 1;
    std::vector<bool> small_sieve(sqrt_max + 1, true);
    small_sieve[0] = small_sieve[1] = false;
    for (int64_t i = 2; i <= sqrt_max; i++)
    {
        if (!small_sieve[i]) continue;
        for (int64_t j = i * i; j <= sqrt_max; j += i)
            small_sieve[j] = false;
    }

    std::vector<bool> is_prime((size_t)range, true);
    for (int64_t i = 2; i <= sqrt_max; i++)
    {
        if (!small_sieve[i]) continue;
        int64_t start = ((min_val + i - 1) / i) * i;
        if (start == i) start += i;
        for (int64_t j = start; j <= max_val; j += i)
            is_prime[(size_t)(j - min_val)] = false;
    }

    for (int64_t i = 0; i < range; i++)
    {
        if (is_prime[(size_t)i])
            result.push_back(min_val + i);
    }
    return result;
}

// ============================================================================
// ABCTemplate: parses and expands ABC-style expression templates
// ============================================================================

struct ABCTemplate
{
    std::vector<std::string> literals;
    std::vector<int> var_indices;
    int num_vars = 0;
    int k_var_index = -1;
    std::string k_fixed;
    std::map<int, int> var_to_pos;
    bool valid = false;

    bool parse(const std::string& header_line)
    {
        valid = false;
        literals.clear();
        var_indices.clear();
        num_vars = 0;
        k_var_index = -1;
        k_fixed.clear();

        if (header_line.size() < 5)
            return false;
        std::string prefix = header_line.substr(0, 4);
        std::string lower_prefix = prefix;
        for (auto& c : lower_prefix) c = (char)std::tolower((unsigned char)c);
        if (lower_prefix != "abc ")
            return false;

        std::string tmpl = header_line.substr(4);
        while (!tmpl.empty() && std::isspace((unsigned char)tmpl.front())) tmpl.erase(tmpl.begin());
        while (!tmpl.empty() && std::isspace((unsigned char)tmpl.back())) tmpl.pop_back();
        if (tmpl.empty())
            return false;

        std::set<int> seen_vars;
        std::string current_literal;
        size_t i = 0;
        while (i < tmpl.size())
        {
            if (tmpl[i] == '$' && i + 1 < tmpl.size() && tmpl[i + 1] >= 'a' && tmpl[i + 1] <= 'd')
            {
                int var_idx = tmpl[i + 1] - 'a';
                seen_vars.insert(var_idx);
                literals.push_back(current_literal);
                current_literal.clear();
                var_indices.push_back(var_idx);
                i += 2;
            }
            else
            {
                current_literal += tmpl[i];
                i++;
            }
        }
        literals.push_back(current_literal);

        if (var_indices.empty())
            return false;

        num_vars = (int)seen_vars.size();
        determine_k_variable(tmpl);

        var_to_pos.clear();
        int pos = 0;
        for (int sv : seen_vars)
            var_to_pos[sv] = pos++;

        valid = true;
        return true;
    }

    bool expand(const std::string& data_line, std::string& expression, std::string& k_value) const
    {
        if (!valid)
            return false;

        std::vector<std::string> values;
        std::istringstream iss(data_line);
        std::string token;
        while (iss >> token)
            values.push_back(token);

        if ((int)values.size() < num_vars)
            return false;

        expression.clear();
        for (size_t i = 0; i < var_indices.size(); i++)
        {
            expression += literals[i];
            expression += get_value(var_indices[i], values);
        }
        expression += literals.back();

        if (k_var_index >= 0)
            k_value = get_value(k_var_index, values);
        else
            k_value = k_fixed;

        return true;
    }

private:
    std::string get_value(int var_idx, const std::vector<std::string>& values) const
    {
        auto it = var_to_pos.find(var_idx);
        if (it != var_to_pos.end() && it->second < (int)values.size())
            return values[it->second];
        return "";
    }

    void determine_k_variable(const std::string& tmpl)
    {
        size_t star_pos = tmpl.find('*');
        if (star_pos == std::string::npos)
        {
            k_var_index = -1;
            k_fixed = "1";
            return;
        }

        std::string k_part = tmpl.substr(0, star_pos);
        for (size_t j = 0; j < k_part.size(); j++)
        {
            if (k_part[j] == '$' && j + 1 < k_part.size() && k_part[j + 1] >= 'a' && k_part[j + 1] <= 'd')
            {
                k_var_index = k_part[j + 1] - 'a';
                return;
            }
        }

        k_var_index = -1;
        k_fixed = k_part;
    }
};

// ============================================================================
// ABC2 variable range parsing
// ============================================================================

struct ABC2VarRange
{
    int var_index;
    std::vector<int64_t> values;
};

static bool parse_abc2_var_line(const std::string& line, ABC2VarRange& vr, Logging& logging, int line_num)
{
    std::string trimmed = str_trim(line);
    if (trimmed.size() < 3 || trimmed[1] != ':')
        return false;

    char var_letter = (char)std::tolower((unsigned char)trimmed[0]);
    if (var_letter < 'a' || var_letter > 'd')
        return false;

    vr.var_index = var_letter - 'a';
    vr.values.clear();

    std::string def = str_trim(trimmed.substr(2));

    // Explicit list: "in { val1 val2 ... }"
    if (def.find("in") != std::string::npos && def.find('{') != std::string::npos)
    {
        size_t brace_start = def.find('{');
        size_t brace_end = def.find('}');
        if (brace_start == std::string::npos || brace_end == std::string::npos || brace_end <= brace_start)
        {
            logging.warning("ABC2: malformed 'in { }' at line %d: %s\n", line_num, line.data());
            return false;
        }
        std::string list_content = def.substr(brace_start + 1, brace_end - brace_start - 1);
        std::istringstream iss(list_content);
        int64_t val;
        while (iss >> val)
            vr.values.push_back(val);
        if (vr.values.size() > MAX_ABC2_IN_LIST_VALUES)
        {
            logging.warning("ABC2: 'in { }' list exceeds %d values at line %d, truncating.\n",
                (int)MAX_ABC2_IN_LIST_VALUES, line_num);
            vr.values.resize(MAX_ABC2_IN_LIST_VALUES);
        }
        return !vr.values.empty();
    }

    // Range definition
    bool use_primes = false;
    bool is_downto = false;
    {
        std::string lower_def = def;
        for (auto& ch : lower_def) ch = (char)std::tolower((unsigned char)ch);
        use_primes = (lower_def.find("primes") != std::string::npos);
        is_downto = (lower_def.find("downto") != std::string::npos);
    }

    // Replace keywords with spaces to extract numbers
    std::string numbers_str = def;
    auto replace_keyword = [&](const char* kw) {
        std::string lower = numbers_str;
        for (auto& ch : lower) ch = (char)std::tolower((unsigned char)ch);
        size_t kwlen = std::strlen(kw);
        std::string lower_kw = kw;
        for (auto& ch : lower_kw) ch = (char)std::tolower((unsigned char)ch);
        size_t pos;
        while ((pos = lower.find(lower_kw)) != std::string::npos)
        {
            for (size_t j = pos; j < pos + kwlen && j < numbers_str.size(); j++)
                numbers_str[j] = ' ';
            for (size_t j = pos; j < pos + kwlen && j < lower.size(); j++)
                lower[j] = ' ';
        }
    };
    replace_keyword("primes");
    replace_keyword("from");
    replace_keyword("downto");
    replace_keyword("to");
    replace_keyword("step");

    std::istringstream iss(numbers_str);
    std::vector<int64_t> nums;
    int64_t num;
    while (iss >> num)
        nums.push_back(num);

    if (nums.size() < 2)
    {
        logging.warning("ABC2: malformed range at line %d: %s\n", line_num, line.data());
        return false;
    }

    int64_t range_start = nums[0];
    int64_t range_end = nums[1];
    int64_t step = (nums.size() >= 3) ? nums[2] : (is_downto ? -1 : 1);

    if (is_downto && step > 0)
    {
        logging.warning("ABC2: 'downto' requires negative step at line %d: %s\n", line_num, line.data());
        return false;
    }
    if (!is_downto && step <= 0)
    {
        logging.warning("ABC2: 'from...to' requires positive step at line %d: %s\n", line_num, line.data());
        return false;
    }
    if (use_primes && nums.size() >= 3)
    {
        logging.warning("ABC2: 'primes' and 'step' cannot be combined at line %d: %s\n", line_num, line.data());
        return false;
    }

    if (use_primes)
    {
        int64_t lo = std::min(range_start, range_end);
        int64_t hi = std::max(range_start, range_end);
        if (hi - lo + 1 > MAX_SIEVE_RANGE)
        {
            logging.warning("ABC2: prime sieve range too large (%lld to %lld, max %lld) at line %d.\n",
                (long long)lo, (long long)hi, (long long)MAX_SIEVE_RANGE, line_num);
            return false;
        }
        vr.values = sieve_primes(lo, hi);
        if (is_downto)
            std::reverse(vr.values.begin(), vr.values.end());
    }
    else if (is_downto)
    {
        for (int64_t v = range_start; v >= range_end; v += step)
            vr.values.push_back(v);
    }
    else
    {
        for (int64_t v = range_start; v <= range_end; v += step)
            vr.values.push_back(v);
    }

    return !vr.values.empty();
}

// ============================================================================
// VectorCandidateSource: materialized vector-backed source
// Used for ABC, ABCD, and raw (non-ABC) formats.
// ============================================================================

class VectorCandidateSource : public CandidateSource
{
public:
    VectorCandidateSource(std::vector<std::string>&& expressions,
                          std::vector<std::string>&& k_values,
                          bool abc)
        : _expressions(std::move(expressions))
        , _k_values(std::move(k_values))
        , _is_abc(abc)
    {
    }

    size_t size() const override { return _expressions.size(); }

    bool get(size_t index, Candidate& out) const override
    {
        if (index >= _expressions.size())
            return false;
        out.expression = _expressions[index];
        out.k_value = (index < _k_values.size()) ? _k_values[index] : "";
        return true;
    }

    bool is_abc() const override { return _is_abc; }

private:
    std::vector<std::string> _expressions;
    std::vector<std::string> _k_values;
    bool _is_abc;
};

// ============================================================================
// ABC2CandidateSource: lazy Cartesian product generator
// Computes candidates on-demand via mixed-radix index decomposition.
// No candidates are stored in memory.
// ============================================================================

class ABC2CandidateSource : public CandidateSource
{
public:
    ABC2CandidateSource(ABCTemplate&& tmpl,
                        std::vector<ABC2VarRange>&& var_ranges,
                        size_t total,
                        std::set<int>&& used_vars,
                        int max_var)
        : _tmpl(std::move(tmpl))
        , _var_ranges(std::move(var_ranges))
        , _total(total)
        , _used_vars(std::move(used_vars))
        , _max_var(max_var)
    {
        // Precompute strides for mixed-radix decomposition
        // Last variable is innermost (fastest changing)
        _strides.resize(_var_ranges.size());
        size_t stride = 1;
        for (int i = (int)_var_ranges.size() - 1; i >= 0; i--)
        {
            _strides[i] = stride;
            stride *= _var_ranges[i].values.size();
        }
    }

    size_t size() const override { return _total; }

    bool get(size_t index, Candidate& out) const override
    {
        if (index >= _total)
            return false;

        // Decompose flat index into per-variable indices
        std::vector<std::string> var_values(_max_var);
        for (size_t r = 0; r < _var_ranges.size(); r++)
        {
            size_t vi = (index / _strides[r]) % _var_ranges[r].values.size();
            var_values[_var_ranges[r].var_index] = std::to_string(_var_ranges[r].values[vi]);
        }

        // Build data line in order expected by ABCTemplate
        std::string data_line;
        for (int vi : _used_vars)
        {
            if (!data_line.empty()) data_line += " ";
            if (vi < (int)var_values.size() && !var_values[vi].empty())
                data_line += var_values[vi];
            else
                data_line += "0";
        }

        return _tmpl.expand(data_line, out.expression, out.k_value);
    }

    bool is_abc() const override { return true; }

private:
    ABCTemplate _tmpl;
    std::vector<ABC2VarRange> _var_ranges;
    size_t _total;
    std::set<int> _used_vars;
    int _max_var;
    std::vector<size_t> _strides;
};

// ============================================================================
// ABCD format parser (materializes into vectors)
// ============================================================================

static bool parse_abcd_file(
    const std::vector<std::string>& raw_lines,
    std::vector<std::string>& out_batch,
    std::vector<std::string>& out_batch_k,
    Logging& logging)
{
    bool any_header = false;
    int header_count = 0;
    ABCTemplate current_template;
    int64_t accum[4] = {0, 0, 0, 0};

    for (size_t i = 0; i < raw_lines.size(); i++)
    {
        const std::string& raw_line = raw_lines[i];
        if (is_skip_line(raw_line))
            continue;

        if (has_prefix_ci(raw_line, "ABCD "))
        {
            std::string header = strip_comment(raw_line);
            size_t bracket_start = header.find('[');
            size_t bracket_end = header.find(']');
            if (bracket_start == std::string::npos || bracket_end == std::string::npos || bracket_end <= bracket_start)
            {
                logging.warning("ABCD: malformed header (missing [...]) at line %d: %s\n", (int)i + 1, raw_line.data());
                continue;
            }

            std::string expr_part = str_trim(header.substr(5, bracket_start - 5));
            std::string abc_header = "ABC " + expr_part;
            if (!current_template.parse(abc_header))
            {
                logging.warning("ABCD: failed to parse expression at line %d: %s\n", (int)i + 1, raw_line.data());
                continue;
            }

            std::string bracket_content = header.substr(bracket_start + 1, bracket_end - bracket_start - 1);
            std::istringstream iss(bracket_content);
            std::string token;
            int var_idx = 0;
            accum[0] = accum[1] = accum[2] = accum[3] = 0;
            while (iss >> token && var_idx < 4)
            {
                try { accum[var_idx] = std::stoll(token); }
                catch (const std::exception&) {
                    logging.warning("ABCD: non-numeric initial value '%s' at line %d\n", token.data(), (int)i + 1);
                    break;
                }
                var_idx++;
            }

            std::string data_line;
            for (int v = 0; v < current_template.num_vars; v++)
            {
                if (v > 0) data_line += " ";
                data_line += std::to_string(accum[v]);
            }
            std::string expression, k_value;
            if (current_template.expand(data_line, expression, k_value))
            {
                out_batch.push_back(std::move(expression));
                out_batch_k.push_back(std::move(k_value));
            }

            any_header = true;
            header_count++;
            continue;
        }

        if (!any_header)
            continue;

        std::istringstream iss(raw_line);
        std::string token;
        int var_idx = 0;
        bool parse_ok = true;
        while (iss >> token && var_idx < current_template.num_vars)
        {
            try {
                int64_t delta = std::stoll(token);
                if ((delta > 0 && accum[var_idx] > INT64_MAX - delta) ||
                    (delta < 0 && accum[var_idx] < INT64_MIN - delta))
                {
                    logging.warning("ABCD: accumulator overflow at line %d, skipping.\n", (int)i + 1);
                    parse_ok = false;
                    break;
                }
                accum[var_idx] += delta;
            }
            catch (const std::exception&) {
                logging.warning("ABCD: non-numeric delta '%s' at line %d, skipping.\n", token.data(), (int)i + 1);
                parse_ok = false;
                break;
            }
            var_idx++;
        }
        if (!parse_ok)
            continue;

        std::string data_line;
        for (int v = 0; v < current_template.num_vars; v++)
        {
            if (v > 0) data_line += " ";
            data_line += std::to_string(accum[v]);
        }
        std::string expression, k_value;
        if (current_template.expand(data_line, expression, k_value))
        {
            out_batch.push_back(std::move(expression));
            out_batch_k.push_back(std::move(k_value));
        }
    }

    if (any_header)
        logging.info("ABCD: parsed %d header blocks, %d candidates.\n", header_count, (int)out_batch.size());
    return any_header && !out_batch.empty();
}

// ============================================================================
// ABC2 format parser (returns lazy source)
// ============================================================================

static std::unique_ptr<CandidateSource> parse_abc2_source(
    const std::vector<std::string>& raw_lines,
    Logging& logging)
{
    if (raw_lines.empty())
        return nullptr;

    std::string header = strip_comment(raw_lines[0]);
    std::string expr_part = str_trim(header.substr(5));
    ABCTemplate tmpl;
    std::string abc_header = "ABC " + expr_part;
    if (!tmpl.parse(abc_header))
    {
        logging.warning("ABC2: failed to parse expression: %s\n", raw_lines[0].data());
        return nullptr;
    }

    logging.info("ABC2 format detected: %s\n", raw_lines[0].data());

    std::vector<ABC2VarRange> var_ranges;
    for (size_t i = 1; i < raw_lines.size(); i++)
    {
        if (is_skip_line(raw_lines[i]))
            continue;
        ABC2VarRange vr;
        if (parse_abc2_var_line(raw_lines[i], vr, logging, (int)i + 1))
            var_ranges.push_back(std::move(vr));
        else
            logging.warning("ABC2: skipping unrecognized line %d: %s\n", (int)i + 1, raw_lines[i].data());
    }

    if (var_ranges.empty())
    {
        logging.warning("ABC2: no valid variable ranges found.\n");
        return nullptr;
    }

    for (const auto& vr : var_ranges)
        logging.info("ABC2: var $%c has %d values\n", (char)('a' + vr.var_index), (int)vr.values.size());

    // Compute total with overflow check
    size_t total = 1;
    bool overflow = false;
    for (const auto& vr : var_ranges)
    {
        if (vr.values.size() > 0 && total > MAX_ABC2_CANDIDATES / vr.values.size())
        {
            overflow = true;
            break;
        }
        total *= vr.values.size();
    }
    if (overflow || total > MAX_ABC2_CANDIDATES)
    {
        logging.warning("ABC2: Cartesian product exceeds %d candidates, aborting.\n", (int)MAX_ABC2_CANDIDATES);
        return nullptr;
    }

    int max_var = 0;
    for (const auto& vr : var_ranges)
        max_var = std::max(max_var, vr.var_index + 1);

    std::set<int> used_vars;
    for (int vi : tmpl.var_indices)
        used_vars.insert(vi);

    logging.info("ABC2: %d candidates (lazy generation).\n", (int)total);

    return std::unique_ptr<CandidateSource>(new ABC2CandidateSource(
        std::move(tmpl), std::move(var_ranges), total, std::move(used_vars), max_var));
}

// ============================================================================
// Factory: parse_batch_file
// ============================================================================

std::unique_ptr<CandidateSource> parse_batch_file(
    const std::string& filename,
    Logging& logging)
{
    File batch_file(filename, 0);
    batch_file.read_buffer();
    std::unique_ptr<TextReader> reader(batch_file.get_textreader());
    std::vector<std::string> raw_lines;
    std::string st;
    while (reader->read_textline(st))
        raw_lines.push_back(std::move(st));

    if (raw_lines.empty())
    {
        logging.warning("Batch %s is empty.\n", filename.data());
        return nullptr;
    }

    FileFormat format = detect_format(raw_lines[0]);

    if (format == FORMAT_ABCD)
    {
        std::vector<std::string> batch, batch_k;
        if (!parse_abcd_file(raw_lines, batch, batch_k, logging))
        {
            logging.warning("ABCD batch %s has no valid data.\n", filename.data());
            return nullptr;
        }
        return std::unique_ptr<CandidateSource>(
            new VectorCandidateSource(std::move(batch), std::move(batch_k), true));
    }

    if (format == FORMAT_ABC2)
    {
        auto source = parse_abc2_source(raw_lines, logging);
        if (!source)
            logging.warning("ABC2 batch %s has no valid data.\n", filename.data());
        return source;
    }

    if (format == FORMAT_ABC)
    {
        ABCTemplate abc_template;
        if (abc_template.parse(raw_lines[0]))
        {
            logging.info("ABC format detected: %s\n", raw_lines[0].data());

            std::vector<std::string> batch, batch_k;
            for (size_t i = 1; i < raw_lines.size(); i++)
            {
                if (is_skip_line(raw_lines[i]))
                    continue;
                std::string expression, k_value;
                if (abc_template.expand(raw_lines[i], expression, k_value))
                {
                    batch.push_back(std::move(expression));
                    batch_k.push_back(std::move(k_value));
                }
                else
                    logging.warning("ABC: skipping malformed data line %d: %s\n", (int)i + 1, raw_lines[i].data());
            }

            if (batch.empty())
            {
                logging.warning("ABC batch %s has no valid data lines.\n", filename.data());
                return nullptr;
            }
            logging.info("ABC: %d candidates from template.\n", (int)batch.size());
            return std::unique_ptr<CandidateSource>(
                new VectorCandidateSource(std::move(batch), std::move(batch_k), true));
        }
    }

    // Raw format: each line is an expression, no k-values
    std::vector<std::string> empty_k;
    return std::unique_ptr<CandidateSource>(
        new VectorCandidateSource(std::move(raw_lines), std::move(empty_k), false));
}
