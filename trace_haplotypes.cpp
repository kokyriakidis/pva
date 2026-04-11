// trace_haplotypes.cpp
// Usage: trace_haplotypes <graph.gbz> <start_node> <end_node>
//
// For each distinct allele spanning the snarl (start_node → end_node),
// extracts the sequence and labels it with all haplotypes that carry it.
//
// Algorithm: iterative DFS on GBWT search states.
//   - No recursion → no stack overflow
//   - Parent pointers for path reconstruction → no string copies during search
//   - index.extend() directly → no CachedGBWT reentrancy issues
//   - follow_edges() for graph topology (collect-first, no callbacks during iteration)

#include <gbwtgraph/gbz.h>
#include <gbwt/gbwt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <unordered_set>
#include <stdexcept>

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
    const gbwt::Metadata& meta = index.metadata;
    gbwt::PathName pname = meta.path(path_id);
    std::string sample = meta.hasSampleNames()
        ? meta.sample(pname.sample)
        : "sample" + std::to_string(pname.sample);
    std::string contig = meta.hasContigNames()
        ? meta.contig(pname.contig)
        : "contig" + std::to_string(pname.contig);
    return sample + "#" + contig
         + "#" + std::to_string(pname.phase)
         + "#" + std::to_string(pname.count);
}

struct Allele {
    std::string       sequence;
    gbwt::SearchState state;
    bool              reversed;
};

// Each visited node in the DFS tree.
struct DFSNode {
    gbwt::SearchState state;
    size_t            parent;   // index in nodes[], SIZE_MAX = root
};

static std::vector<Allele> search(
    const gbwtgraph::GBWTGraph& graph,
    const gbwt::GBWT&           index,
    gbwt::SearchState           start_state,
    gbwt::node_type             end_fwd,
    gbwt::node_type             end_rev,
    bool                        reversed)
{
    std::vector<Allele>  alleles;
    std::vector<DFSNode> nodes;   // DFS tree, grows monotonically
    std::vector<size_t>  stk;     // indices into nodes[]

    nodes.reserve(4096);
    nodes.push_back({start_state, SIZE_MAX});
    stk.push_back(0);

    while (!stk.empty()) {
        size_t idx = stk.back(); stk.pop_back();

        // Copy state — avoid dangling ref if nodes[] reallocates below
        gbwt::SearchState state = nodes[idx].state;

        // Collect graph successors (non-recursive — no reentrancy issue)
        handlegraph::handle_t h = gbwtgraph::GBWTGraph::node_to_handle(state.node);
        std::vector<handlegraph::handle_t> next_handles;
        graph.follow_edges(h, false, [&](handlegraph::handle_t nh) {
            next_handles.push_back(nh);
            return true;
        });

        for (handlegraph::handle_t nh : next_handles) {
            gbwt::node_type next_node = gbwt::Node::encode(
                graph.get_id(nh), graph.get_is_reverse(nh));

            // Extend GBWT state — only haplotypes that actually cross this edge survive
            gbwt::SearchState next_state = index.extend(state, next_node);
            if (next_state.empty()) continue;

            size_t child_idx = nodes.size();
            nodes.push_back({next_state, idx});

            if (next_node == end_fwd || next_node == end_rev) {
                // Reconstruct path by walking parent chain
                std::vector<gbwt::node_type> path_nodes;
                for (size_t i = child_idx; i != SIZE_MAX; i = nodes[i].parent)
                    path_nodes.push_back(nodes[i].state.node);
                std::reverse(path_nodes.begin(), path_nodes.end());

                std::string seq;
                seq.reserve(path_nodes.size() * 32);
                for (gbwt::node_type n : path_nodes)
                    seq += graph.get_sequence(gbwtgraph::GBWTGraph::node_to_handle(n));

                alleles.push_back({std::move(seq), next_state, reversed});
                // Don't push to stack — end node reached
            } else {
                stk.push_back(child_idx);
            }
        }
    }
    return alleles;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <graph.gbz> <start_node> <end_node>\n";
        return 1;
    }

    gbwtgraph::nid_t start_id = std::stoll(argv[2]);
    gbwtgraph::nid_t end_id   = std::stoll(argv[3]);

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

    for (auto [search_node, reversed] :
         std::vector<std::pair<gbwt::node_type, bool>>{
             {start_fwd, false}, {start_rev, true}}) {

        if (!index.contains(search_node)) continue;
        gbwt::SearchState state = index.find(search_node);
        if (state.empty()) continue;

        auto found = search(graph, index, state, end_fwd, end_rev, reversed);
        alleles.insert(alleles.end(),
                       std::make_move_iterator(found.begin()),
                       std::make_move_iterator(found.end()));
    }

    std::unordered_set<gbwt::size_type> seen;
    for (const auto& allele : alleles) {
        std::vector<gbwt::size_type> seq_ids = index.locate(allele.state);
        for (gbwt::size_type seq_id : seq_ids) {
            gbwt::size_type path_id = gbwt::Path::id(seq_id);
            if (!seen.insert(path_id).second) continue;

            std::string seq = allele.reversed
                ? reverse_complement(allele.sequence)
                : allele.sequence;

            std::cout << ">" << get_label(index, seq_id) << "\n"
                      << seq << "\n";
        }
    }

    return 0;
}
