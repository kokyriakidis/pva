// pva.h — Pangenome VNTR Aligner public API
//
// Link against libpva.so to call these functions from your own C++ program,
// or use the `pva` binary for command-line access.
//
// Each function returns 0 on success, non-zero on failure.
//
// Note: pva::trace() requires linking against gbwtgraph + gbwt + sdsl.
//       All other functions are in libpva.so with no heavy dependencies.

#pragma once

#include <string>
#include <cstdint>
#include <vector>

namespace pva {

// ── query ─────────────────────────────────────────────────────────────────────
// Extract a subgraph GFA for a genomic interval from a GBZ database.
// Wraps: gbz-base/target/release/query
int query(const std::string& db,
          const std::string& sample,
          const std::string& contig,
          uint64_t           interval_start,
          uint64_t           interval_end,
          const std::string& out_gfa,
          const std::string& query_bin = "");

// ── snarls ────────────────────────────────────────────────────────────────────
// Find snarls in a GFA subgraph.
// Wraps: vg snarls <gfa> > <out_snarls>
int snarls(const std::string& gfa,
           const std::string& out_snarls,
           const std::string& vg_bin = "");

// ── boundary ──────────────────────────────────────────────────────────────────
// Find the genomically outermost boundary nodes of all top-level snarls.
// Uses the reference sample's W line from the GFA to determine genomic order.
struct BoundaryNodes {
    int64_t start_node = -1;
    int64_t end_node   = -1;
};
int boundary(const std::string& gfa,
             const std::string& snarls_file,
             const std::string& sample,
             BoundaryNodes&     out,
             const std::string& vg_bin = "");

// ── guide_tree ────────────────────────────────────────────────────────────────
// Prepare all-vs-all centrolign run for guide tree inference via SLURM.
// Splits multi-FASTA, generates combinations.txt, and writes SLURM commands.
// After SLURM jobs finish, run `pva infer-guide-tree` to get the Newick tree.
struct SlurmOptions {
    std::string partition;
    std::string mail_user;
    std::string mail_type;
    std::string mem;
    std::string time_limit;
    std::string centrolign_bin;
    int         nodes          = 1;
    int         cpus_per_task  = 0;
    std::string log_pattern    = "logs/array_job_%A_task_%a.log";
    std::string job_name_prefix = "vntr_all_pairs";
};
int guide_tree(const std::string& fasta,
               const std::string& outdir,
               const std::string& slurm_out,
               const SlurmOptions& slurm,
               int                max_array  = 30000);

// ── align ─────────────────────────────────────────────────────────────────────
// Pairwise alignment via centrolign (in-memory, no subprocess).
// Returns explicit CIGAR string (=/X ops).
// align_pair / align_all_pairs are in libpva.so — requires libcentrolign.so.
std::string pairwise_cigar(const std::string& name1, const std::string& seq1,
                            const std::string& name2, const std::string& seq2,
                            bool skip_calibration = false);

int align_pair(const std::string& name1, const std::string& seq1,
               const std::string& name2, const std::string& seq2,
               const std::string& out_path,
               bool skip_calibration = false);

int align_all_pairs(const std::string& fasta,
                    const std::string& cigar_dir,
                    bool skip_calibration = false);

// ── infer_tree ────────────────────────────────────────────────────────────────
// Infer a neighbor-joining guide tree from pairwise CIGAR files.
// Wraps scripts/infer_tree.py via fork+execvp.
int infer_guide_tree(const std::string& cigar_dir,
               const std::string& combinations,
               const std::string& out_path,
               const std::string& script_hint = "");

// ── snarl_dist_filter ─────────────────────────────────────────────────────────
// A boundary node pair extracted from a snarl.
struct SnarlRegion {
    int64_t start_node = -1;
    int64_t end_node   = -1;
};

// Parse all top-level snarls from a binary snarls file and return their
// boundary node pairs. Uses `vg view -R` internally.
std::vector<SnarlRegion> read_snarl_regions(const std::string& snarls_file,
                                             const std::string& vg_bin = "");

// Filter top-level snarls by minimum distance between their boundary nodes,
// queried from the distance index built by `vg index -j`.
// Output: filtered binary snarls file (re-encoded via `vg view -JR -`).
int dist_filter(const std::string& snarls_file,
                const std::string& dist_file,
                const std::string& out_path,
                int64_t            min_dist    = 0,
                int64_t            max_dist    = INT64_MAX,
                const std::string& vg_bin      = "");

// ── CLI entry points ──────────────────────────────────────────────────────────
// Each *_main() receives its own argc/argv (argv[0] == subcommand name).
// trace_main() is NOT in libpva.so — it requires gbwtgraph linkage.
// align_main() requires libcentrolign.so linkage.
int query_main       (int argc, char* argv[]);
int snarls_main      (int argc, char* argv[]);
int boundary_main    (int argc, char* argv[]);
int prep_guide_tree_slurm_main (int argc, char* argv[]);
int snarl_dist_filter_main (int argc, char* argv[]);
int infer_guide_tree_main        (int argc, char* argv[]);
int trace_main             (int argc, char* argv[]); // linked into pva binary only
int align_main             (int argc, char* argv[]); // linked into pva binary only

} // namespace pva
