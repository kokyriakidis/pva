// pva_distfilter.cpp — `pva snarl-dist-filter` subcommand
//
// Filters top-level snarls by the minimum distance between their boundary nodes,
// as reported by `vg distance` using the snarl distance index.
//
// For each top-level snarl in the snarls file:
//   1. Extract its start and end boundary node IDs from `vg view -R` JSON
//   2. Query `vg distance -d <dist_file> -s <start>:0 -e <end>:0`
//   3. Keep snarls whose minimum distance falls within [--min-dist, --max-dist]
//
// Output: a filtered binary snarls file, re-encoded by piping the passing
// NDJSON lines through `vg view -JR -`. This file can be fed directly to
// `pva trace --snarls` or any other tool that consumes vg snarls.

#include "pva.h"
#include "pva_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <climits>

namespace pva {

// ── Snarl parsing ─────────────────────────────────────────────────────────────

static bool line_is_top_level(const std::string& line) {
    auto pos = line.find("\"parent\"");
    if (pos == std::string::npos) return true;
    auto colon = line.find(':', pos);
    if (colon == std::string::npos) return true;
    size_t v = colon + 1;
    while (v < line.size() && line[v] == ' ') ++v;
    return line.substr(v, 4) == "null";
}

// Extract start.node_id and end.node_id from a snarl JSON line.
static bool parse_snarl_region(const std::string& line, SnarlRegion& out) {
    static const std::regex re_start(
        R"re("start"\s*:\s*\{[^}]*"node_id"\s*:\s*"?(\d+)"?)re");
    static const std::regex re_end(
        R"re("end"\s*:\s*\{[^}]*"node_id"\s*:\s*"?(\d+)"?)re");
    std::smatch sm, em;
    if (!std::regex_search(line, sm, re_start)) return false;
    if (!std::regex_search(line, em, re_end))   return false;
    out.start_node = std::stoll(sm[1].str());
    out.end_node   = std::stoll(em[1].str());
    return true;
}

// Internal: a parsed snarl entry keeping the original JSON line for re-encoding.
struct SnarlEntry {
    SnarlRegion region;
    std::string json_line;
};

static std::vector<SnarlEntry> parse_snarl_entries(const std::string& ndjson) {
    std::vector<SnarlEntry> entries;
    std::istringstream ss(ndjson);
    std::string line;
    while (std::getline(ss, line)) {
        if (line.empty() || !line_is_top_level(line)) continue;
        SnarlEntry e;
        if (parse_snarl_region(line, e.region)) {
            e.json_line = line;
            entries.push_back(std::move(e));
        }
    }
    return entries;
}

// ── Distance query ────────────────────────────────────────────────────────────

static const int64_t DIST_UNREACHABLE = INT64_MAX;
static const int64_t VG_MAX_DIST      = (int64_t)9223372036854775807LL;

static int64_t query_distance(const std::string& vg_bin,
                               const std::string& dist_file,
                               int64_t start_node, int64_t end_node)
{
    std::string sp = std::to_string(start_node) + ":0";
    std::string ep = std::to_string(end_node)   + ":0";

    std::string output;
    int rc = pva::utils::run_capture(
        {vg_bin, "distance", "-d", dist_file, "-s", sp, "-e", ep}, output);
    if (rc != 0) return -1;

    // Extract the last integer from output (vg prints "... : 42" or just "42")
    std::smatch m;
    std::regex re_int(R"(\b(\d+)\b)");
    auto it = output.cbegin();
    int64_t dist = -1;
    while (std::regex_search(it, output.cend(), m, re_int)) {
        dist = std::stoll(m[1].str());
        it = m.suffix().first;
    }
    if (dist < 0)          return -1;
    if (dist >= VG_MAX_DIST) return DIST_UNREACHABLE;
    return dist;
}

// ── Public API ────────────────────────────────────────────────────────────────

std::vector<SnarlRegion> read_snarl_regions(const std::string& snarls_file,
                                             const std::string& vg_hint)
{
    std::string vg_bin = pva::utils::find_bin(vg_hint, "vg", "vg");
    std::string ndjson;
    if (pva::utils::run_capture({vg_bin, "view", "-R", snarls_file}, ndjson) != 0) {
        std::cerr << "pva: vg view -R failed for " << snarls_file << "\n";
        return {};
    }
    auto entries = parse_snarl_entries(ndjson);
    std::vector<SnarlRegion> regions;
    regions.reserve(entries.size());
    for (auto& e : entries) regions.push_back(e.region);
    return regions;
}

int dist_filter(const std::string& snarls_file,
                const std::string& dist_file,
                const std::string& out_path,
                int64_t            min_d,
                int64_t            max_d,
                const std::string& vg_hint)
{
    std::string vg_bin = pva::utils::find_bin(vg_hint, "vg", "vg");

    std::string ndjson;
    if (pva::utils::run_capture({vg_bin, "view", "-R", snarls_file}, ndjson) != 0) {
        std::cerr << "pva snarl-dist-filter: vg view -R failed\n"; return 1;
    }

    auto entries = parse_snarl_entries(ndjson);
    if (entries.empty()) {
        std::cerr << "pva snarl-dist-filter: no top-level snarls parsed\n"; return 1;
    }

    // Filter by distance, collecting passing JSON lines for re-encoding
    std::string filtered_json;
    size_t kept = 0, skipped_unreachable = 0, errors = 0;
    int ret = 0;

    for (const auto& e : entries) {
        int64_t d = query_distance(vg_bin, dist_file,
                                   e.region.start_node, e.region.end_node);
        if (d < 0) {
            std::cerr << "pva snarl-dist-filter: distance query failed for snarl "
                      << e.region.start_node << " -> " << e.region.end_node << "\n";
            ++errors; ret = 1; continue;
        }
        if (d == DIST_UNREACHABLE) { ++skipped_unreachable; continue; }
        if (d >= min_d && d <= max_d) {
            filtered_json += e.json_line + "\n";
            ++kept;
        }
    }

    std::string max_str = (max_d == INT64_MAX) ? "inf" : std::to_string(max_d);
    std::cerr << "pva snarl-dist-filter: kept " << kept << "/" << entries.size()
              << " snarls  [min=" << min_d << " max=" << max_str << "]";
    if (skipped_unreachable > 0)
        std::cerr << "  (" << skipped_unreachable << " unreachable)";
    if (errors > 0)
        std::cerr << "  (" << errors << " errors)";
    std::cerr << "\n";

    if (kept == 0) {
        std::cerr << "pva snarl-dist-filter: no snarls passed filter\n";
        return ret;
    }

    // Re-encode filtered NDJSON → binary snarls via `vg view -JR -`
    int rc = pva::utils::run_pipe_to_file(
        {vg_bin, "view", "-JR", "-"}, filtered_json, out_path);
    if (rc != 0) {
        std::cerr << "pva snarl-dist-filter: vg view -JR failed (exit " << rc << ")\n";
        return 1;
    }
    return ret;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void distfilter_usage(const char* prog) {
    std::cerr <<
        "Usage:\n"
        "  " << prog << " snarl-dist-filter [options]\n"
        "\n"
        "Filter top-level snarls by the minimum distance between their boundary nodes.\n"
        "Uses `vg distance` with the snarl distance index built by `vg index -j`.\n"
        "Output is a filtered binary snarls file consumable by `pva trace --snarls`.\n"
        "\n"
        "Required:\n"
        "  --snarls <path>   Input snarls file\n"
        "  --dist <path>     Snarl distance index file (.dist)\n"
        "  --out <path>      Output filtered snarls file\n"
        "\n"
        "Optional:\n"
        "  --min-dist <N>    Keep snarls with min distance >= N (default: 0)\n"
        "  --max-dist <N>    Keep snarls with min distance <= N (default: unlimited)\n"
        "\n"
        "Unreachable snarl pairs (disconnected nodes) are excluded from output.\n"
        "\n"
        "Example:\n"
        "  pva snarl-dist-filter --snarls region.snarls --dist graph.dist \\\n"
        "                        --out filtered.snarls --min-dist 100 --max-dist 50000\n"
        "  pva trace --gbz graph.gbz --snarls filtered.snarls --out-dir alleles/\n";
}

int snarl_dist_filter_main(int argc, char* argv[]) {
    int64_t min_d  = 0;
    int64_t max_d  = INT64_MAX;
    std::string snarls_file, dist_file, out_path, vg_hint;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva snarl-dist-filter: " << a << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--snarls")            snarls_file = next();
        else if (a == "--dist")              dist_file   = next();
        else if (a == "--out")               out_path    = next();
        else if (a == "--min-dist")          min_d       = std::stoll(next());
        else if (a == "--max-dist")          max_d       = std::stoll(next());
        else if (a == "--vg-bin")            vg_hint     = next();
        else if (a == "--help" || a == "-h") { distfilter_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva snarl-dist-filter: unknown option " << a << "\n";
            distfilter_usage(argv[0]); return 1;
        }
    }

    if (snarls_file.empty() || dist_file.empty() || out_path.empty()) {
        distfilter_usage(argv[0]); return 1;
    }

    return pva::dist_filter(snarls_file, dist_file, out_path, min_d, max_d, vg_hint);
}

} // namespace pva
