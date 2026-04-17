// pva_trace_haplotypes.cpp — `pva trace-haplotypes` subcommand
//
// Traces all haplotype sequences spanning a snarl (start_node → end_node)
// using the GBWT r-index for fast locate().
//
// This is the integrated haplotype-tracing implementation used by the
// `pva trace-haplotypes` subcommand. Requires gbwtgraph + gbwt + sdsl.

#include "pva.h"

#include <gbwtgraph/gbz.h>
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <cerrno>
#include <sys/stat.h>
#include <linux/limits.h>

namespace pva {

// ── helpers ───────────────────────────────────────────────────────────────────

static std::string reverse_complement(const std::string& s) {
    std::string rc(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            default:  rc[i] = c;   break;
        }
    }
    return rc;
}

static std::string get_label(const gbwt::GBWT& index, gbwt::size_type seq_id) {
    gbwt::size_type path_id = gbwt::Path::id(seq_id);
    if (!index.hasMetadata() || !index.metadata.hasPathNames())
        return "path_" + std::to_string(path_id);
    gbwt::FullPathName name = index.metadata.fullPath(path_id);
    return name.sample_name + "#" + name.contig_name
         + "#" + std::to_string(name.haplotype)
         + "#" + std::to_string(name.offset);
}

// ── r-index load/build ────────────────────────────────────────────────────────

static gbwt::FastLocate build_or_load_rindex(const gbwt::GBWT& index,
                                              const std::string& gbz_path,
                                              const std::string& ri_hint = "")
{
    std::string ri_path = ri_hint.empty() ? gbz_path + ".ri" : ri_hint;
    {
        std::ifstream f(ri_path, std::ios::binary);
        if (f) {
            std::cerr << "Loading r-index from " << ri_path << "...\n";
            gbwt::FastLocate r;
            r.load(f);
            r.setGBWT(index);
            std::cerr << "r-index loaded.\n";
            return r;
        }
    }
    std::cerr << "Building r-index (will be cached as " << ri_path << ")...\n";
    gbwt::FastLocate r(index);
    {
        std::ofstream f(ri_path, std::ios::binary);
        if (f) r.serialize(f);
        else   std::cerr << "Warning: could not cache r-index to " << ri_path << "\n";
    }
    std::cerr << "r-index built.\n";
    return r;
}

// ── DFS ───────────────────────────────────────────────────────────────────────

struct DFSNode {
    gbwt::SearchState state;
    gbwt::size_type   first;
    size_t            parent;
};

struct Allele {
    std::string       sequence;
    gbwt::SearchState state;
    gbwt::size_type   first;
};

static std::vector<Allele> dfs_search(
    const gbwtgraph::GBWTGraph& graph,
    const gbwt::FastLocate&     r_index,
    gbwt::SearchState           start_state,
    gbwt::size_type             start_first,
    gbwt::node_type             end_fwd,
    gbwt::node_type             end_rev)
{
    std::vector<Allele>  alleles;
    std::vector<DFSNode> nodes;
    std::vector<size_t>  stk;

    nodes.reserve(4096);
    nodes.push_back({start_state, start_first, SIZE_MAX});
    stk.push_back(0);

    while (!stk.empty()) {
        size_t idx = stk.back(); stk.pop_back();
        gbwt::SearchState state = nodes[idx].state;
        gbwt::size_type   first = nodes[idx].first;

        handlegraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(state.node);
        std::vector<handlegraph::handle_t> next_handles;
        graph.follow_edges(h, false, [&](handlegraph::handle_t nh) {
            next_handles.push_back(nh); return true;
        });

        for (handlegraph::handle_t nh : next_handles) {
            gbwt::node_type next_node = gbwt::Node::encode(
                graph.get_id(nh), graph.get_is_reverse(nh));

            gbwt::size_type next_first = first;
            gbwt::SearchState next_state = r_index.extend(state, next_node, next_first);
            if (next_state.empty()) continue;

            size_t child_idx = nodes.size();
            nodes.push_back({next_state, next_first, idx});

            if (next_node == end_fwd || next_node == end_rev) {
                std::vector<gbwt::node_type> path_nodes;
                for (size_t i = child_idx; i != SIZE_MAX; i = nodes[i].parent)
                    path_nodes.push_back(nodes[i].state.node);
                std::reverse(path_nodes.begin(), path_nodes.end());

                std::string seq;
                seq.reserve(path_nodes.size() * 32);
                for (gbwt::node_type n : path_nodes)
                    seq += graph.get_sequence(gbwtgraph::GBWTGraph::node_to_handle(n));

                alleles.push_back({std::move(seq), next_state, next_first});
            } else {
                stk.push_back(child_idx);
            }
        }
    }
    return alleles;
}

// ── Single region ─────────────────────────────────────────────────────────────

static int trace_region(
    const gbwtgraph::GBWTGraph& graph,
    const gbwt::GBWT&           index,
    const gbwt::FastLocate&     r_index,
    gbwtgraph::nid_t            start_id,
    gbwtgraph::nid_t            end_id,
    std::ostream&               out)
{
    if (!graph.has_node(start_id)) {
        std::cerr << "pva trace-haplotypes: start node " << start_id << " not in graph\n"; return 1;
    }
    if (!graph.has_node(end_id)) {
        std::cerr << "pva trace-haplotypes: end node " << end_id << " not in graph\n"; return 1;
    }

    gbwt::node_type start_fwd = gbwt::Node::encode(start_id, false);
    gbwt::node_type start_rev = gbwt::Node::encode(start_id, true);
    gbwt::node_type end_fwd   = gbwt::Node::encode(end_id,   false);
    gbwt::node_type end_rev   = gbwt::Node::encode(end_id,   true);

    std::vector<Allele> alleles;
    for (gbwt::node_type search_node : {start_fwd, start_rev}) {
        if (!index.contains(search_node)) continue;
        gbwt::size_type first = gbwt::FastLocate::NO_POSITION;
        gbwt::SearchState state = r_index.find(search_node, first);
        if (state.empty()) continue;
        auto found = dfs_search(graph, r_index, state, first, end_fwd, end_rev);
        alleles.insert(alleles.end(),
                       std::make_move_iterator(found.begin()),
                       std::make_move_iterator(found.end()));
    }

    std::unordered_set<gbwt::size_type> seen_paths;
    for (const auto& allele : alleles) {
        std::vector<gbwt::size_type> seq_ids = r_index.locate(allele.state, allele.first);
        for (gbwt::size_type seq_id : seq_ids) {
            gbwt::size_type path_id = gbwt::Path::id(seq_id);
            if (!seen_paths.insert(path_id).second) continue;
            std::string seq = gbwt::Path::is_reverse(seq_id)
                ? reverse_complement(allele.sequence)
                : allele.sequence;
            out << ">" << get_label(index, seq_id) << "\n" << seq << "\n";
        }
    }
    return 0;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void trace_usage(const char* prog) {
    std::cerr <<
        "Usage:\n"
        "  " << prog << " trace --gbz <graph.gbz> --start <node> --end <node> [--out <out.fa>]\n"
        "  " << prog << " trace --gbz <graph.gbz> --snarls <filtered.snarls> --out-dir <dir>\n"
        "\n"
        "Options:\n"
        "  --gbz <path>      Input GBZ graph file\n"
        "  --start <node>    Start boundary node ID (single-region mode)\n"
        "  --end <node>      End boundary node ID (single-region mode)\n"
        "  --out <path>      Output FASTA (single-region mode, default: stdout)\n"
        "  --snarls <path>   Filtered snarls file (output of pva snarl-dist-filter)\n"
        "  --out-dir <path>  Output directory for per-snarl FASTAs (required with --snarls)\n"
        "  --ri <path>       Pre-built r-index file (default: <graph.gbz>.ri, built and cached on first run)\n"
        "\n"
        "The r-index (~15-50x faster locate()) is built automatically on first run\n"
        "and cached as <graph.gbz>.ri. Use --ri to specify a custom path.\n"
        "\n"
        "Example:\n"
        "  pva trace-haplotypes --gbz graph.gbz --start 12345 --end 12360 --out alleles.fa\n"
        "  pva trace-haplotypes --gbz graph.gbz --snarls filtered.snarls --out-dir alleles/\n";
}

int trace_main(int argc, char* argv[]) {
    std::string gbz_path, snarls_file, out_path, out_dir, ri_path;
    gbwtgraph::nid_t start_id = 0, end_id = 0;
    bool has_start = false, has_end = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva trace-haplotypes: " << a << " requires an argument\n"; std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--gbz")               gbz_path    = next();
        else if (a == "--start")           { start_id    = std::stoll(next()); has_start = true; }
        else if (a == "--end")             { end_id      = std::stoll(next()); has_end   = true; }
        else if (a == "--out")               out_path    = next();
        else if (a == "--snarls")            snarls_file = next();
        else if (a == "--out-dir")           out_dir     = next();
        else if (a == "--ri")                ri_path     = next();
        else if (a == "--help" || a == "-h") { trace_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva trace-haplotypes: unknown option " << a << "\n";
            trace_usage(argv[0]); return 1;
        }
    }

    if (gbz_path.empty()) { trace_usage(argv[0]); return 1; }

    // Load GBZ
    gbwtgraph::GBZ gbz;
    try {
        std::ifstream in(gbz_path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open " + gbz_path);
        gbz.simple_sds_load(in);
    } catch (const std::exception& e) {
        std::cerr << "pva trace-haplotypes: " << e.what() << "\n"; return 1;
    }
    const gbwt::GBWT&           index = gbz.index;
    const gbwtgraph::GBWTGraph& graph = gbz.graph;

    gbwt::FastLocate r_index = build_or_load_rindex(index, gbz_path, ri_path);

    // Snarls mode: read filtered snarls, auto-name output FASTAs
    if (!snarls_file.empty()) {
        if (out_dir.empty()) {
            std::cerr << "pva trace-haplotypes: --snarls requires --out-dir\n";
            trace_usage(argv[0]); return 1;
        }
        // Ensure out_dir exists
        struct stat st;
        if (stat(out_dir.c_str(), &st) != 0) {
            if (mkdir(out_dir.c_str(), 0755) != 0 && errno != EEXIST) {
                std::cerr << "pva trace-haplotypes: cannot create " << out_dir << "\n"; return 1;
            }
        }

        auto regions = pva::read_snarl_regions(snarls_file);
        if (regions.empty()) {
            std::cerr << "pva trace-haplotypes: no regions in " << snarls_file << "\n"; return 1;
        }

        int ret = 0;
        size_t done = 0;
        for (const auto& r : regions) {
            std::string opath = out_dir + "/snarl_"
                              + std::to_string(r.start_node) + "_"
                              + std::to_string(r.end_node)   + ".fa";
            std::cerr << "[" << ++done << "/" << regions.size() << "] "
                      << r.start_node << " -> " << r.end_node
                      << " => " << opath << "\n";
            std::ofstream out_file(opath);
            if (!out_file) {
                std::cerr << "pva trace-haplotypes: cannot open " << opath << "\n";
                ret = 1; continue;
            }
            ret |= trace_region(graph, index, r_index,
                                (gbwtgraph::nid_t)r.start_node,
                                (gbwtgraph::nid_t)r.end_node, out_file);
        }
        return ret;
    }

    // Single mode
    if (!has_start || !has_end) { trace_usage(argv[0]); return 1; }

    if (!out_path.empty() && out_path != "-") {
        std::ofstream out_file(out_path);
        if (!out_file) {
            std::cerr << "pva trace-haplotypes: cannot open " << out_path << "\n"; return 1;
        }
        return trace_region(graph, index, r_index, start_id, end_id, out_file);
    }
    return trace_region(graph, index, r_index, start_id, end_id, std::cout);
}

} // namespace pva
