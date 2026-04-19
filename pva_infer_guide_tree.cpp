// pva_infer_tree.cpp — `pva infer-guide-tree` subcommand
//
// Infers a neighbor-joining guide tree from pairwise centrolign CIGAR files.
// Wraps scripts/infer_tree.py via fork+execvp — the NJ algorithm is implemented
// there using scikit-bio (install via: pip install scikit-bio).
//
// The script is auto-detected relative to the pva binary:
//   <pva_dir>/../scripts/infer_tree.py
//
// Inputs come from the output directory of `pva prep-guide-tree-slurm`:
//   --cigar-dir <path>    Directory containing pairwise_cigar_*.txt files
//                         produced by the completed SLURM jobs
//   --combinations <path> combinations.txt created during prep
//   --out <path>          Output Newick file (default: stdout)

#include "pva.h"
#include "pva_utils.h"

#include <iostream>
#include <string>
#include <vector>
#include <linux/limits.h>
#include <unistd.h>

namespace pva {

static std::string bin_dir() {
    char self[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (len <= 0) return {};
    self[len] = '\0';
    std::string path(self);
    auto slash = path.rfind('/');
    return slash != std::string::npos ? path.substr(0, slash) : std::string{};
}

static std::string find_python() {
    std::string dir = bin_dir();
    if (!dir.empty()) {
        std::string venv_py = dir + "/../venv/bin/python3";
        if (access(venv_py.c_str(), X_OK) == 0)
            return venv_py;
    }
    return "python3";
}

static std::string find_infer_tree_script(const std::string& hint) {
    if (!hint.empty()) return hint;

    std::string dir = bin_dir();
    if (!dir.empty()) {
        std::string candidate = dir + "/../scripts/infer_tree.py";
        if (access(candidate.c_str(), R_OK) == 0)
            return candidate;
    }
    return "scripts/infer_tree.py";
}

// ── Public API ────────────────────────────────────────────────────────────────

int infer_guide_tree(const std::string& cigar_dir,
               const std::string& combinations,
               const std::string& out_path,
               const std::string& script_hint)
{
    std::string script = find_infer_tree_script(script_hint);
    std::vector<std::string> args = {find_python(), script, cigar_dir, combinations};

    if (out_path.empty() || out_path == "-") {
        // Let stdout flow through to the terminal / caller's stdout
        std::vector<const char*> argv;
        for (auto& s : args) argv.push_back(s.c_str());
        argv.push_back(nullptr);

        pid_t pid = fork();
        if (pid < 0) {
            std::cerr << "pva infer-guide-tree: fork failed\n"; return 1;
        }
        if (pid == 0) {
            execvp(argv[0], const_cast<char* const*>(argv.data()));
            std::cerr << "pva infer-guide-tree: exec failed\n"; _exit(127);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        if (!WIFEXITED(status)) return 1;
        return WEXITSTATUS(status);
    }

    return pva::utils::run_to_file(args, out_path);
}

// ── CLI ───────────────────────────────────────────────────────────────────────

static void infer_tree_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " infer-guide-tree [options]\n"
        "\n"
        "Infer a neighbor-joining guide tree from pairwise centrolign CIGAR files.\n"
        "Inputs are the output of `pva prep-guide-tree-slurm` + completed SLURM alignment jobs.\n"
        "\n"
        "Required:\n"
        "  --cigar-dir <path>     Directory containing pairwise_cigar_*.txt files from completed SLURM jobs\n"
        "  --combinations <path>  combinations.txt created by pva prep-guide-tree-slurm\n"
        "  --out <path>           Output Newick file\n"
        "\n"
        "Optional:\n"
        "  --script <path>        Path to infer_tree.py (auto-detected if omitted)\n"
        "\n"
        "Example:\n"
        "  pva infer-guide-tree \\\n"
        "      --cigar-dir guide_tree/pairwise_cigar \\\n"
        "      --combinations guide_tree/combinations.txt \\\n"
        "      --out guide_tree/guide_tree.nwk\n";
}

int infer_guide_tree_main(int argc, char* argv[]) {
    std::string cigar_dir, combinations, out_path, script_hint;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva infer-guide-tree: " << a << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if      (a == "--cigar-dir")         cigar_dir    = next();
        else if (a == "--combinations")      combinations = next();
        else if (a == "--out")               out_path     = next();
        else if (a == "--script")            script_hint  = next();
        else if (a == "--help" || a == "-h") { infer_tree_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva infer-guide-tree: unknown option " << a << "\n";
            infer_tree_usage(argv[0]); return 1;
        }
    }

    if (cigar_dir.empty() || combinations.empty() || out_path.empty()) {
        infer_tree_usage(argv[0]); return 1;
    }

    return pva::infer_guide_tree(cigar_dir, combinations, out_path, script_hint);
}

} // namespace pva
