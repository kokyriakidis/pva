// pva_boundary.cpp — `pva boundary` subcommand
//
// Finds the genomically outermost boundary nodes of all top-level snarls in a
// GFA subgraph, using the reference sample's walk to determine genomic order
// (node IDs are not in genomic order in pangenome graphs).
//
// Steps:
//   1. Run `vg view -R <snarls>` and capture JSON output
//   2. Parse each JSON line: collect start/end node_ids for top-level snarls
//      (those without a "parent" field, or with "parent":null)
//   3. Scan the GFA for the first W line belonging to <sample>, walk it left
//      to right; first boundary node seen = genomic start, last = genomic end

#include "pva.h"
#include "pva_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <unordered_set>
#include <regex>
#include <cctype>

namespace pva {

// ── JSON parsing helpers ──────────────────────────────────────────────────────
// vg view -R outputs one protobuf-to-JSON Snarl object per line (NDJSON).
// We only need to check two things per line:
//   - Is it a top-level snarl? (no "parent":{} field, or "parent":null)
//   - What are the start and end node_ids?

static bool is_top_level(const std::string& line) {
    auto pos = line.find("\"parent\"");
    if (pos == std::string::npos) return true;   // field absent → top-level
    // field present — top-level only if value is null
    auto colon = line.find(':', pos);
    if (colon == std::string::npos) return true;
    size_t v = colon + 1;
    while (v < line.size() && line[v] == ' ') ++v;
    return line.substr(v, 4) == "null";
}

static std::unordered_set<int64_t> extract_node_ids(const std::string& line) {
    // Matches "node_id":"12345" or "node_id":12345
    static const std::regex re(R"(\"node_id\"\s*:\s*\"?(\d+)\"?)");
    std::unordered_set<int64_t> ids;
    auto begin = std::sregex_iterator(line.begin(), line.end(), re);
    auto end   = std::sregex_iterator();
    for (auto it = begin; it != end; ++it)
        ids.insert(std::stoll((*it)[1].str()));
    return ids;
}

// ── GFA W-line walker ─────────────────────────────────────────────────────────
// Finds the first W line for <sample> and walks it in order, tracking which
// boundary nodes appear first and last.
static bool find_boundary_in_gfa(const std::string&                 gfa_path,
                                  const std::string&                 sample,
                                  const std::unordered_set<int64_t>& boundary,
                                  int64_t&                           start_out,
                                  int64_t&                           end_out,
                                  std::string&                       error_out)
{
    std::ifstream gfa(gfa_path);
    if (!gfa) {
        error_out = "pva boundary: cannot open GFA: " + gfa_path;
        return false;
    }

    std::string line;
    bool saw_w_line = false;
    std::set<std::string> sample_names;
    while (std::getline(gfa, line)) {
        // W lines: W  <sample>  <hap>  <contig>  <start>  <end>  <walk>
        if (line.size() < 2 || line[0] != 'W' || line[1] != '\t') continue;
        saw_w_line = true;

        // Check sample name (field 2)
        size_t f1 = line.find('\t');                       // after 'W'
        size_t f2 = line.find('\t', f1 + 1);              // after sample
        if (f1 == std::string::npos || f2 == std::string::npos) continue;
        std::string line_sample = line.substr(f1 + 1, f2 - f1 - 1);
        sample_names.insert(line_sample);
        if (line_sample != sample) continue;

        // Walk field is the last tab-separated field
        size_t walk_start = line.rfind('\t');
        if (walk_start == std::string::npos) continue;
        std::string walk = line.substr(walk_start + 1);

        // Parse walk: e.g. ">12344>12345<12347>12350"
        size_t pos = 0;
        while (pos < walk.size()) {
            if (walk[pos] == '>' || walk[pos] == '<') {
                size_t num_start = pos + 1;
                size_t num_end   = num_start;
                while (num_end < walk.size() && std::isdigit((unsigned char)walk[num_end]))
                    ++num_end;
                if (num_end > num_start) {
                    int64_t node = std::stoll(walk.substr(num_start, num_end - num_start));
                    if (boundary.count(node)) {
                        if (start_out == -1) start_out = node;
                        end_out = node;
                    }
                }
                pos = num_end;
            } else {
                ++pos;
            }
        }
        return true;  // only process the first matching W line
    }

    if (!saw_w_line) {
        error_out = "pva boundary: no W lines found in " + gfa_path;
        return false;
    }

    std::ostringstream msg;
    msg << "pva boundary: sample '" << sample << "' not found in any W line of "
        << gfa_path;
    if (!sample_names.empty()) {
        msg << " (available samples:";
        bool first = true;
        for (const auto& name : sample_names) {
            msg << (first ? " " : ", ") << name;
            first = false;
        }
        msg << ")";
    }
    error_out = msg.str();
    return false;
}

// ── Public API ────────────────────────────────────────────────────────────────

int boundary(const std::string& gfa,
             const std::string& snarls_file,
             const std::string& sample,
             BoundaryNodes&     out,
             const std::string& vg_bin)
{
    out.start_node = -1;
    out.end_node   = -1;

    // 1. Run vg view -R to get JSON snarls
    std::string bin = pva::utils::find_bin(vg_bin, "vg", "vg");
    std::string json_output;
    int ret = pva::utils::run_capture({bin, "view", "-R", snarls_file}, json_output);
    if (ret != 0) {
        std::cerr << "pva boundary: vg view -R failed (exit " << ret << ")\n";
        return 1;
    }

    // 2. Parse JSON lines to collect top-level boundary node IDs
    std::unordered_set<int64_t> boundary_nodes;
    std::istringstream ss(json_output);
    std::string line;
    while (std::getline(ss, line)) {
        if (line.empty()) continue;
        if (!is_top_level(line)) continue;
        auto ids = extract_node_ids(line);
        boundary_nodes.insert(ids.begin(), ids.end());
    }

    if (boundary_nodes.empty()) {
        std::cerr << "pva boundary: no top-level snarl boundary nodes found\n";
        return 1;
    }

    // 3. Find genomic order from the GFA W line
    std::string gfa_error;
    if (!find_boundary_in_gfa(gfa, sample, boundary_nodes,
                               out.start_node, out.end_node, gfa_error)) {
        std::cerr << gfa_error << "\n";
        return 1;
    }

    if (out.start_node == -1 || out.end_node == -1) {
        std::cerr << "pva boundary: boundary nodes not found in " << sample
                  << " walk\n";
        return 1;
    }
    if (out.start_node == out.end_node) {
        std::cerr << "pva boundary: start and end node are the same ("
                  << out.start_node << ") — region may be too small\n";
        return 1;
    }

    return 0;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void boundary_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " boundary [options]\n"
        "\n"
        "Options:\n"
        "  --gfa <path>      Input GFA file\n"
        "  --snarls <path>   Input snarls file\n"
        "  --sample <name>   Reference sample used to determine genomic order (required; must exist in a GFA W line)\n"
        "  --vg-bin <path>   Path to vg binary (auto-detected if omitted)\n"
        "\n"
        "Output:\n"
        "  start_node end_node   (tab-separated, written to stdout)\n"
        "\n"
        "Example:\n"
        "  pva boundary --gfa region.gfa --snarls region.snarls --sample CHM13\n";
}

int boundary_main(int argc, char* argv[]) {
    std::string gfa, snarls_file, sample, vg_bin;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva boundary: " << a << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--gfa")               gfa         = next();
        else if (a == "--snarls")            snarls_file = next();
        else if (a == "--sample")            sample      = next();
        else if (a == "--vg-bin")            vg_bin      = next();
        else if (a == "--help" || a == "-h") { boundary_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva boundary: unknown option " << a << "\n";
            boundary_usage(argv[0]); return 1;
        }
    }

    if (gfa.empty() || snarls_file.empty() || sample.empty()) {
        if (sample.empty()) {
            std::cerr << "pva boundary: --sample is required\n";
        }
        boundary_usage(argv[0]);
        return 1;
    }

    BoundaryNodes nodes;
    int ret = pva::boundary(gfa, snarls_file, sample, nodes, vg_bin);
    if (ret != 0) return ret;

    std::cout << nodes.start_node << "\t" << nodes.end_node << "\n";
    return 0;
}

} // namespace pva
