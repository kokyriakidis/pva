// pva_align.cpp — `pva align` subcommand
//
// Runs centrolign pairwise alignment using libcentrolign.so directly.
// No subprocess — sequences are passed in-memory, CIGAR returned as string.
//
// For a pair of sequences:
//   1. Build a trivial two-leaf guide tree via centrolign::in_order_newick_string()
//   2. Construct centrolign::Core in-memory (no temp files)
//   3. core.execute()
//   4. Induce pairwise alignment from the root subproblem graph
//   5. Emit explicit CIGAR (=/X ops) to stdout or file
//
// For all-vs-all from a multi-FASTA:
//   - Parse FASTA once into memory
//   - Generate N*(N-1)/2 pairs, run each pair sequentially
//   - Write CIGAR to pairwise_cigar/<name1>_<name2>.txt  (mirrors guide_tree pipeline)

#include "pva.h"

#include <centrolign/core.hpp>
#include <centrolign/tree.hpp>
#include <centrolign/alignment.hpp>
#include <centrolign/logging.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include <sys/stat.h>

namespace pva {

// ── helpers ───────────────────────────────────────────────────────────────────

static std::string safe_name(const std::string& label) {
    return std::regex_replace(label, std::regex("[^0-9a-zA-Z.\\-]+"), "_");
}

static bool mkdir_p(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) return true;
    auto slash = path.rfind('/');
    if (slash != std::string::npos && slash > 0)
        if (!mkdir_p(path.substr(0, slash))) return false;
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
}

struct FastaEntry { std::string label; std::string sequence; };

static std::vector<FastaEntry> read_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "pva align: cannot open " << path << "\n"; return {}; }
    std::vector<FastaEntry> entries;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') entries.push_back({line.substr(1), {}});
        else if (!entries.empty()) entries.back().sequence += line;
    }
    return entries;
}

// ── core pairwise alignment ───────────────────────────────────────────────────

std::string pairwise_cigar(const std::string& name1, const std::string& seq1,
                            const std::string& name2, const std::string& seq2,
                            bool skip_calibration)
{
    // Trivial two-leaf Newick: (name1,name2);
    std::string newick = centrolign::in_order_newick_string({name1, name2});
    centrolign::Tree tree(newick);

    std::vector<std::pair<std::string, std::string>> seqs = {
        {name1, seq1},
        {name2, seq2}
    };

    centrolign::Core core(std::move(seqs), std::move(tree));
    core.skip_calibration = skip_calibration;

    core.execute();

    // Induce pairwise alignment from the completed root subproblem graph
    const auto& graph = core.root_subproblem().graph;

    uint64_t path_id1 = graph.path_id(name1);
    uint64_t path_id2 = graph.path_id(name2);

    centrolign::Alignment aln = centrolign::induced_pairwise_alignment(graph, path_id1, path_id2);

    return centrolign::explicit_cigar(aln, graph, graph);
}

// ── public API ────────────────────────────────────────────────────────────────

int align_pair(const std::string& name1, const std::string& seq1,
               const std::string& name2, const std::string& seq2,
               const std::string& out_path,
               bool skip_calibration)
{
    std::string cigar = pairwise_cigar(name1, seq1, name2, seq2, skip_calibration);

    if (out_path.empty() || out_path == "-") {
        std::cout << cigar << "\n";
    } else {
        std::ofstream f(out_path);
        if (!f) {
            std::cerr << "pva align: cannot write " << out_path << "\n";
            return 1;
        }
        f << cigar << "\n";
    }
    return 0;
}

int align_all_pairs(const std::string& fasta,
                    const std::string& cigar_dir,
                    bool skip_calibration)
{
    if (!mkdir_p(cigar_dir)) {
        std::cerr << "pva align: cannot create " << cigar_dir << "\n";
        return 1;
    }

    auto entries = read_fasta(fasta);
    if (entries.size() < 2) {
        std::cerr << "pva align: need at least 2 sequences\n";
        return 1;
    }

    size_t n = entries.size();
    size_t n_pairs = n * (n - 1) / 2;
    size_t done = 0;

    int ret = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const auto& e1 = entries[i];
            const auto& e2 = entries[j];

            std::string s1 = safe_name(e1.label);
            std::string s2 = safe_name(e2.label);
            // lexicographic order to match combinations.txt convention
            if (s1 > s2) std::swap(s1, s2);

            std::string out = cigar_dir + "/pairwise_cigar_" + s1 + "_" + s2 + ".txt";

            // Skip already-completed pairs (resumable)
            {
                struct stat st;
                if (stat(out.c_str(), &st) == 0 && st.st_size > 0) {
                    ++done;
                    continue;
                }
            }

            std::cerr << "[" << ++done << "/" << n_pairs << "] "
                      << e1.label << " vs " << e2.label << "\n";

            ret |= align_pair(e1.label, e1.sequence, e2.label, e2.sequence,
                              out, skip_calibration);
        }
    }
    return ret;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void align_usage(const char* prog) {
    std::cerr <<
        "Usage:\n"
        "  " << prog << " align --seq1 <seq1.fa> --seq2 <seq2.fa> [--out <cigar>]\n"
        "  " << prog << " align --all-pairs --fasta <seqs.fa> --out <cigar_dir>\n"
        "\n"
        "Options:\n"
        "  --seq1 <path>        First FASTA file (single-pair mode)\n"
        "  --seq2 <path>        Second FASTA file (single-pair mode)\n"
        "  --all-pairs          Run all-vs-all pairwise alignment\n"
        "  --fasta <path>       Multi-FASTA input (all-pairs mode)\n"
        "  --out <path>         Output CIGAR file or directory (default: stdout)\n"
        "  --skip-calibration   Skip anchor score calibration (faster, use for guide tree)\n"
        "  -v <level>           Verbosity: 0=silent 1=minimal 2=basic 3=verbose (default: 0)\n"
        "\n"
        "Output:\n"
        "  Explicit CIGAR string (=/X ops) to stdout or --out file.\n"
        "  --all-pairs writes pairwise_cigar_<name1>_<name2>.txt into <--out dir>/.\n"
        "\n"
        "Example:\n"
        "  pva align --seq1 HG002.fa --seq2 HG003.fa\n"
        "  pva align --all-pairs --fasta alleles.fa --out cigar_dir/ --skip-calibration\n";
}

int align_main(int argc, char* argv[]) {
    bool all_pairs        = false;
    bool skip_calibration = false;
    std::string seq1, seq2, fasta, out_path;
    int verbosity = 0;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva align: " << a << " requires an argument\n"; std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--seq1")              seq1             = next();
        else if (a == "--seq2")              seq2             = next();
        else if (a == "--fasta")             fasta            = next();
        else if (a == "--out")               out_path         = next();
        else if (a == "--all-pairs")         all_pairs        = true;
        else if (a == "--skip-calibration")  skip_calibration = true;
        else if (a == "-v")                  verbosity        = std::stoi(next());
        else if (a == "--help" || a == "-h") { align_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva align: unknown option " << a << "\n";
            align_usage(argv[0]); return 1;
        }
    }

    // Set centrolign verbosity
    switch (verbosity) {
        case 0: centrolign::logging::level = centrolign::logging::Silent;  break;
        case 1: centrolign::logging::level = centrolign::logging::Minimal; break;
        case 2: centrolign::logging::level = centrolign::logging::Basic;   break;
        default: centrolign::logging::level = centrolign::logging::Verbose; break;
    }

    if (all_pairs) {
        if (fasta.empty() || out_path.empty()) { align_usage(argv[0]); return 1; }
        return pva::align_all_pairs(fasta, out_path, skip_calibration);
    }

    if (seq1.empty() || seq2.empty()) { align_usage(argv[0]); return 1; }

    auto e1 = read_fasta(seq1);
    auto e2 = read_fasta(seq2);
    if (e1.empty() || e2.empty()) {
        std::cerr << "pva align: could not read sequences\n"; return 1;
    }

    return pva::align_pair(e1[0].label, e1[0].sequence,
                           e2[0].label, e2[0].sequence,
                           out_path, skip_calibration);
}

} // namespace pva
