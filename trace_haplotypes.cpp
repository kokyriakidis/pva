// trace_haplotypes.cpp
//
// Usage (single region):
//   trace_haplotypes <graph.gbz> <start_node> <end_node>
//
// Usage (batch — load GBZ once, query many regions):
//   trace_haplotypes <graph.gbz> --batch <regions.tsv>
//
// regions.tsv format (tab-separated, one region per line):
//   start_node  end_node  output.fa
//
// Algorithm (per region):
//   - Iterative DFS on GBWT BidirectionalState (bdFind + bdExtendForward)
//   - bdExtendForward: O(1) amortized — narrows haplotype set to those crossing each edge
//   - Dead branches (no haplotype takes that edge) are pruned immediately
//   - Parent pointers for path reconstruction — no string copies during search
//   - index.locate() maps each terminal state to individual haplotype IDs

#include <gbwtgraph/gbz.h>
#include <gbwt/gbwt.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>

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
    // fullPath() is the official API: resolves sample/contig indices to strings
    // and returns a FullPathName with all fields populated in one call.
    gbwt::FullPathName name = index.metadata.fullPath(path_id);
    return name.sample_name + "#" + name.contig_name
         + "#" + std::to_string(name.haplotype)
         + "#" + std::to_string(name.offset);
}

// ── DFS types ─────────────────────────────────────────────────────────────────

struct DFSNode {
    gbwt::BidirectionalState state;
    size_t                   parent;  // index in nodes[], SIZE_MAX = root
};

struct Allele {
    std::string              sequence;
    gbwt::BidirectionalState state;    // at end_node — for locate()
    // No 'reversed' field: orientation is read from seq_id via Path::is_reverse()
    // after locate(), which is the authoritative source in the GBWT metadata.
};

// ── core search ───────────────────────────────────────────────────────────────

static std::vector<Allele> search(
    const gbwtgraph::GBWTGraph& graph,
    const gbwt::GBWT&           index,
    gbwt::BidirectionalState    start_state,
    gbwt::node_type             end_fwd,
    gbwt::node_type             end_rev)
{
    std::vector<Allele>  alleles;
    std::vector<DFSNode> nodes;
    std::vector<size_t>  stk;

    nodes.reserve(4096);
    nodes.push_back({start_state, SIZE_MAX});
    stk.push_back(0);

    while (!stk.empty()) {
        size_t idx = stk.back(); stk.pop_back();

        gbwt::BidirectionalState state = nodes[idx].state;

        handlegraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(state.forward.node);
        std::vector<handlegraph::handle_t> next_handles;
        graph.follow_edges(h, false, [&](handlegraph::handle_t nh) {
            next_handles.push_back(nh);
            return true;
        });

        for (handlegraph::handle_t nh : next_handles) {
            gbwt::node_type next_node = gbwt::Node::encode(
                graph.get_id(nh), graph.get_is_reverse(nh));

            // bdExtendForward: O(1) amortized; returns empty if no haplotype
            // crosses this edge — dead branches pruned immediately
            gbwt::BidirectionalState next_state = index.bdExtendForward(state, next_node);
            if (next_state.empty()) continue;

            size_t child_idx = nodes.size();
            nodes.push_back({next_state, idx});

            if (next_node == end_fwd || next_node == end_rev) {
                std::vector<gbwt::node_type> path_nodes;
                for (size_t i = child_idx; i != SIZE_MAX; i = nodes[i].parent)
                    path_nodes.push_back(nodes[i].state.forward.node);
                std::reverse(path_nodes.begin(), path_nodes.end());

                std::string seq;
                seq.reserve(path_nodes.size() * 32);
                for (gbwt::node_type n : path_nodes)
                    seq += graph.get_sequence(gbwtgraph::GBWTGraph::node_to_handle(n));

                alleles.push_back({std::move(seq), next_state});
            } else {
                stk.push_back(child_idx);
            }
        }
    }
    return alleles;
}

// ── emit one region ───────────────────────────────────────────────────────────

static int trace_region(
    const gbwtgraph::GBWTGraph& graph,
    const gbwt::GBWT&           index,
    gbwtgraph::nid_t            start_id,
    gbwtgraph::nid_t            end_id,
    std::ostream&               out)
{
    if (!graph.has_node(start_id)) {
        std::cerr << "Error: start node " << start_id << " not in graph\n";
        return 1;
    }
    if (!graph.has_node(end_id)) {
        std::cerr << "Error: end node " << end_id << " not in graph\n";
        return 1;
    }

    gbwt::node_type start_fwd = gbwt::Node::encode(start_id, false);
    gbwt::node_type start_rev = gbwt::Node::encode(start_id, true);
    gbwt::node_type end_fwd   = gbwt::Node::encode(end_id,   false);
    gbwt::node_type end_rev   = gbwt::Node::encode(end_id,   true);

    std::vector<Allele> alleles;
    for (gbwt::node_type search_node : {start_fwd, start_rev}) {
        if (!index.contains(search_node)) continue;
        gbwt::BidirectionalState state = index.bdFind(search_node);
        if (state.empty()) continue;

        auto found = search(graph, index, state, end_fwd, end_rev);
        alleles.insert(alleles.end(),
                       std::make_move_iterator(found.begin()),
                       std::make_move_iterator(found.end()));
    }

    std::unordered_set<gbwt::size_type> seen_paths;
    for (const auto& allele : alleles) {
        std::vector<gbwt::size_type> seq_ids = index.locate(allele.state.forward);
        for (gbwt::size_type seq_id : seq_ids) {
            gbwt::size_type path_id = gbwt::Path::id(seq_id);
            if (!seen_paths.insert(path_id).second) continue;

            // Path::is_reverse() reads orientation directly from the seq_id —
            // the authoritative source, no need to track it through the DFS.
            std::string seq = gbwt::Path::is_reverse(seq_id)
                ? reverse_complement(allele.sequence)
                : allele.sequence;

            out << ">" << get_label(index, seq_id) << "\n" << seq << "\n";
        }
    }
    return 0;
}

// ── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <graph.gbz> <start_node> <end_node>\n"
                  << "  " << argv[0] << " <graph.gbz> --batch <regions.tsv>\n"
                  << "\nregions.tsv: start_node\\tend_node\\toutput.fa  (one per line)\n";
        return 1;
    }

    // ── Load GBZ once ─────────────────────────────────────────────────────────
    gbwtgraph::GBZ gbz;
    try {
        std::ifstream in(argv[1], std::ios::binary);
        if (!in) throw std::runtime_error("cannot open file");
        gbz.simple_sds_load(in);
    } catch (const std::exception& e) {
        std::cerr << "Error loading GBZ: " << e.what() << "\n";
        return 1;
    }
    const gbwt::GBWT&           index = gbz.index;
    const gbwtgraph::GBWTGraph& graph = gbz.graph;

    // ── Batch mode ─────────────────────────────────────────────────────────────
    if (std::string(argv[2]) == "--batch") {
        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " <graph.gbz> --batch <regions.tsv>\n";
            return 1;
        }
        std::ifstream tsv(argv[3]);
        if (!tsv) {
            std::cerr << "Error: cannot open " << argv[3] << "\n";
            return 1;
        }

        std::string line;
        int ret = 0;
        while (std::getline(tsv, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            gbwtgraph::nid_t start_id, end_id;
            std::string out_path;
            if (!(ss >> start_id >> end_id >> out_path)) {
                std::cerr << "Warning: skipping malformed line: " << line << "\n";
                continue;
            }
            std::cerr << "Tracing " << start_id << " -> " << end_id
                      << " => " << out_path << "\n";
            std::ofstream out(out_path);
            if (!out) {
                std::cerr << "Error: cannot open output " << out_path << "\n";
                ret = 1;
                continue;
            }
            ret |= trace_region(graph, index, start_id, end_id, out);
        }
        return ret;
    }

    // ── Single mode ────────────────────────────────────────────────────────────
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <graph.gbz> <start_node> <end_node>\n";
        return 1;
    }
    gbwtgraph::nid_t start_id = std::stoll(argv[2]);
    gbwtgraph::nid_t end_id   = std::stoll(argv[3]);
    return trace_region(graph, index, start_id, end_id, std::cout);
}
