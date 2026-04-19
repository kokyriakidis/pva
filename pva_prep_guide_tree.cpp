// pva_guide_tree.cpp — `pva guide-tree` subcommand
//
// Prepares all-vs-all centrolign run for guide tree inference via SLURM:
//   1. Splits multi-FASTA into per-haplotype FASTA files (split_fastas/)
//   2. Builds fasta_list.txt  (absolute paths, one per line, sorted)
//   3. Generates combinations.txt  (path1 \t path2, lexicographic, N*(N-1)/2 pairs)
//   4. Writes SLURM sbatch commands to a file (chunked to respect --max-array limit)

#include "pva.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <sys/stat.h>
#include <unistd.h>
#include <linux/limits.h>
#include <cerrno>

namespace pva {

// ── helpers ───────────────────────────────────────────────────────────────────

static std::string safe_name(const std::string& label) {
    return std::regex_replace(label, std::regex("[^0-9a-zA-Z.\\-]+"), "_");
}

static bool mkdir_p(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) return true;
    // Create parent first
    auto slash = path.rfind('/');
    if (slash != std::string::npos && slash > 0)
        if (!mkdir_p(path.substr(0, slash))) return false;
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
}

static std::string self_dir() {
    char self[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (len <= 0) return ".";
    self[len] = '\0';
    std::string p(self);
    auto slash = p.rfind('/');
    return (slash != std::string::npos) ? p.substr(0, slash) : ".";
}

// ── FASTA splitting ───────────────────────────────────────────────────────────

struct FastaEntry { std::string label; std::string sequence; };

static std::vector<FastaEntry> read_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "pva guide-tree: cannot open " << path << "\n";
        return {};
    }
    std::vector<FastaEntry> entries;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            entries.push_back({line.substr(1), {}});
        } else if (!entries.empty()) {
            entries.back().sequence += line;
        }
    }
    return entries;
}

static int split_fasta(const std::vector<FastaEntry>& entries,
                        const std::string& split_dir)
{
    for (const auto& e : entries) {
        std::string path = split_dir + "/" + safe_name(e.label) + ".fa";
        std::ofstream out(path);
        if (!out) {
            std::cerr << "pva guide-tree: cannot write " << path << "\n";
            return 1;
        }
        out << ">" << safe_name(e.label) << "\n" << e.sequence << "\n";
    }
    return 0;
}

// ── Public API ────────────────────────────────────────────────────────────────

int guide_tree(const std::string& fasta,
               const std::string& outdir,
               const std::string& slurm_out,
               const SlurmOptions& slurm_cfg,
               int                max_array)
{
    std::string split_dir    = outdir + "/split_fastas";
    std::string cigar_dir    = outdir + "/pairwise_cigar";
    std::string log_dir      = outdir + "/logs";
    std::string fasta_list   = outdir + "/fasta_list.txt";
    std::string combinations = outdir + "/combinations.txt";

    for (auto& d : {outdir, split_dir, cigar_dir, log_dir}) {
        if (!mkdir_p(d)) {
            std::cerr << "pva guide-tree: cannot create directory " << d << "\n";
            return 1;
        }
    }

    // 1. Read and split multi-FASTA
    std::cerr << "Splitting haplotypes...\n";
    auto entries = read_fasta(fasta);
    if (entries.empty()) {
        std::cerr << "pva guide-tree: no sequences in " << fasta << "\n";
        return 1;
    }
    if (split_fasta(entries, split_dir) != 0) return 1;
    std::cerr << "Haplotypes: " << entries.size() << "\n";

    if (entries.size() < 2) {
        std::cerr << "pva guide-tree: need at least 2 haplotypes\n";
        return 1;
    }

    // 2. Build fasta_list.txt (absolute paths, sorted)
    std::vector<std::string> fa_paths;
    fa_paths.reserve(entries.size());
    for (const auto& e : entries) {
        // Resolve absolute path like realpath
        std::string rel = split_dir + "/" + safe_name(e.label) + ".fa";
        char abs[PATH_MAX];
        if (realpath(rel.c_str(), abs))
            fa_paths.push_back(abs);
        else
            fa_paths.push_back(rel);
    }
    std::sort(fa_paths.begin(), fa_paths.end());

    {
        std::ofstream out(fasta_list);
        if (!out) { std::cerr << "pva guide-tree: cannot write " << fasta_list << "\n"; return 1; }
        for (auto& p : fa_paths) out << p << "\n";
    }

    // 3. Generate all-vs-all combinations (lexicographic, no duplicate pairs)
    std::cerr << "Generating combinations...\n";
    size_t n = fa_paths.size();
    size_t n_pairs = n * (n - 1) / 2;

    {
        std::ofstream out(combinations);
        if (!out) { std::cerr << "pva guide-tree: cannot write " << combinations << "\n"; return 1; }
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i + 1; j < n; ++j)
                if (fa_paths[i] < fa_paths[j])
                    out << fa_paths[i] << "\t" << fa_paths[j] << "\n";
                else
                    out << fa_paths[j] << "\t" << fa_paths[i] << "\n";
    }
    std::cerr << "Total pairs: " << n_pairs << "\n";

    // 4. Write SLURM commands (chunked to respect max_array limit)
    std::string bin_dir    = self_dir();
    std::string slurm_script = bin_dir + "/../slurm/all_pairs.sh";

    // Resolve combinations absolute path
    char abs_comb[PATH_MAX];
    std::string combinations_abs = realpath(combinations.c_str(), abs_comb)
                                 ? std::string(abs_comb) : combinations;
    char abs_outdir[PATH_MAX];
    std::string outdir_abs = realpath(outdir.c_str(), abs_outdir)
                           ? std::string(abs_outdir) : outdir;

    auto slash = slurm_out.rfind('/');
    if (slash != std::string::npos && slash > 0) {
        if (!mkdir_p(slurm_out.substr(0, slash))) {
            std::cerr << "pva guide-tree: cannot create directory "
                      << slurm_out.substr(0, slash) << "\n";
            return 1;
        }
    }

    std::ofstream slurm(slurm_out);
    if (!slurm) {
        std::cerr << "pva guide-tree: cannot write " << slurm_out << "\n";
        return 1;
    }

    slurm << "cd " << outdir_abs << "\n\n";

    size_t start = 1;
    int    batch = 1;
    while (start <= n_pairs) {
        size_t end = std::min(start + (size_t)max_array - 1, n_pairs);
        slurm << "# Batch " << batch << " (tasks " << start << "-" << end << ")\n"
              << "sbatch \\\n"
              << "    --partition=" << slurm_cfg.partition << " \\\n"
              << "    --nodes=" << slurm_cfg.nodes << " \\\n"
              << "    --mem=" << slurm_cfg.mem << " \\\n"
              << "    --ntasks=1 \\\n"
              << "    --cpus-per-task=" << slurm_cfg.cpus_per_task << " \\\n"
              << "    --output=" << slurm_cfg.log_pattern << " \\\n"
              << "    --time=" << slurm_cfg.time_limit << " \\\n";
        if (!slurm_cfg.mail_user.empty()) {
            slurm << "    --mail-user=" << slurm_cfg.mail_user << " \\\n";
        }
        if (!slurm_cfg.mail_type.empty()) {
            slurm << "    --mail-type=" << slurm_cfg.mail_type << " \\\n";
        }
        slurm << "    --job-name=" << slurm_cfg.job_name_prefix << "_" << batch << " \\\n"
              << "    --array=[" << start << "-" << end << "]%" << slurm_cfg.cpus_per_task << " \\\n"
              << "    --export=COMBINATIONS_FILE=" << combinations_abs
              << ",OUTDIR=" << outdir_abs
              << ",CENTROLIGN_BIN=" << slurm_cfg.centrolign_bin << " \\\n"
              << "    " << slurm_script << "\n\n";
        start = end + 1;
        ++batch;
    }

    std::cerr << "Wrote SLURM commands to " << slurm_out << "\n";

    return 0;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void guide_tree_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " prep-guide-tree-slurm [options]\n"
        "\n"
        "Options:\n"
        "  --fasta <path>     Input multi-FASTA of haplotype sequences\n"
        "  --out-dir <path>   Output directory for split FASTAs and combinations.txt\n"
        "  --slurm-out <path> Output file that will contain SLURM submission commands\n"
        "  --partition <name> SLURM partition/queue name\n"
        "  --mem <spec>       SLURM memory request (e.g. 56gb)\n"
        "  --time <spec>      SLURM time limit (e.g. 1:00:00)\n"
        "  --centrolign-bin <path> Path to centrolign binary used by the SLURM worker\n"
        "  --cpus-per-task <n> SLURM CPUs per task and array parallelism cap\n"
        "  --nodes <n>        SLURM node count (default: 1)\n"
        "  --log-pattern <path> SLURM log path pattern relative to --out-dir (default: logs/array_job_%A_task_%a.log)\n"
        "  --job-name-prefix <name> Prefix for generated SLURM job names (default: vntr_all_pairs)\n"
        "  --mail-user <addr> Email address for SLURM notifications (optional)\n"
        "  --mail-type <type> SLURM mail type (optional; e.g. ALL, END, FAIL)\n"
        "  --max-array <n>    SLURM array job size limit (default: 30000)\n"
        "\n"
        "Example:\n"
        "  pva prep-guide-tree-slurm --fasta data/alleles.fa --out-dir data/guide_tree --slurm-out data/guide_tree/submit_jobs.sh \\\n"
        "      --partition short --mem 56gb --time 1:00:00 --cpus-per-task 32 \\\n"
        "      --centrolign-bin /path/to/centrolign\n";
}

int prep_guide_tree_slurm_main(int argc, char* argv[]) {
    std::string fasta, outdir, slurm_out;
    SlurmOptions slurm_cfg;
    int max_array  = 30000;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << "pva guide-tree: " << a << " requires an argument\n"; std::exit(1); }
            return argv[++i];
        };
        auto next_int = [&](int& v) { v = std::stoi(next()); };
        if      (a == "--fasta")             fasta  = next();
        else if (a == "--out-dir")           outdir = next();
        else if (a == "--slurm-out")         slurm_out = next();
        else if (a == "--partition")         slurm_cfg.partition = next();
        else if (a == "--mail-user")         slurm_cfg.mail_user = next();
        else if (a == "--mail-type")         slurm_cfg.mail_type = next();
        else if (a == "--mem")               slurm_cfg.mem = next();
        else if (a == "--time")              slurm_cfg.time_limit = next();
        else if (a == "--centrolign-bin")    slurm_cfg.centrolign_bin = next();
        else if (a == "--nodes")             next_int(slurm_cfg.nodes);
        else if (a == "--cpus-per-task")     next_int(slurm_cfg.cpus_per_task);
        else if (a == "--log-pattern")       slurm_cfg.log_pattern = next();
        else if (a == "--job-name-prefix")   slurm_cfg.job_name_prefix = next();
        else if (a == "--max-array")         next_int(max_array);
        else if (a == "--help" || a == "-h") { guide_tree_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva guide-tree: unknown option " << a << "\n";
            guide_tree_usage(argv[0]); return 1;
        }
    }

    if (fasta.empty() || outdir.empty() || slurm_out.empty() ||
        slurm_cfg.partition.empty() || slurm_cfg.mem.empty() ||
        slurm_cfg.time_limit.empty() || slurm_cfg.centrolign_bin.empty() ||
        slurm_cfg.cpus_per_task <= 0) {
        guide_tree_usage(argv[0]);
        return 1;
    }
    return pva::guide_tree(fasta, outdir, slurm_out, slurm_cfg, max_array);
}

} // namespace pva
