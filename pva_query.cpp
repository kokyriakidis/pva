// pva_query.cpp — implementation of the `pva query` subcommand
//
// Wraps gbz-base/target/release/query. The gbz-base query binary is the
// authoritative tool for extracting GFA subgraphs from .db files; we call it
// via fork+execvp with a stdout redirect so argument quoting is exact and
// there is no shell injection surface.

#include "pva.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cerrno>
#include <climits>

#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <linux/limits.h>

namespace pva {

// ── binary discovery ──────────────────────────────────────────────────────────
// Returns the path to the gbz-base query binary.
// Preference order:
//   1. Caller-supplied path
//   2. <dir_of_pva_binary>/gbz-base/target/release/query
//   3. "query" on PATH
static std::string find_query_bin(const std::string& hint) {
    if (!hint.empty()) return hint;

    // Find directory of the currently running pva binary via /proc/self/exe
    char self[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (len > 0) {
        self[len] = '\0';
        std::string dir(self);
        auto slash = dir.rfind('/');
        if (slash != std::string::npos) {
            std::string candidate = dir.substr(0, slash + 1)
                                  + "gbz-base/target/release/query";
            if (access(candidate.c_str(), X_OK) == 0)
                return candidate;
        }
    }

    return "query"; // fall back to PATH
}

// ── core implementation ───────────────────────────────────────────────────────
int query(const std::string& db,
          const std::string& sample,
          const std::string& contig,
          uint64_t           interval_start,
          uint64_t           interval_end,
          const std::string& out_gfa,
          const std::string& query_bin)
{
    std::string bin = find_query_bin(query_bin);
    std::string interval = std::to_string(interval_start)
                         + ".."
                         + std::to_string(interval_end);

    // Build argv — no shell, no quoting issues
    std::vector<std::string> args = {
        bin,
        "--sample",   sample,
        "--contig",   contig,
        "--interval", interval,
        db
    };
    std::vector<const char*> argv;
    for (auto& s : args) argv.push_back(s.c_str());
    argv.push_back(nullptr);

    // Open output file before fork so we can report errors in the parent
    int fd = open(out_gfa.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        std::cerr << "pva query: cannot open output " << out_gfa
                  << ": " << strerror(errno) << "\n";
        return 1;
    }

    pid_t pid = fork();
    if (pid < 0) {
        std::cerr << "pva query: fork failed: " << strerror(errno) << "\n";
        close(fd);
        return 1;
    }

    if (pid == 0) {
        // child: redirect stdout → output file, then exec
        if (dup2(fd, STDOUT_FILENO) < 0) _exit(1);
        close(fd);
        execvp(argv[0], const_cast<char* const*>(argv.data()));
        // execvp only returns on failure
        std::cerr << "pva query: exec failed: " << strerror(errno)
                  << " (tried: " << bin << ")\n";
        _exit(127);
    }

    close(fd);

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        std::cerr << "pva query: waitpid failed: " << strerror(errno) << "\n";
        return 1;
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::cerr << "pva query: gbz-base query exited with status "
                  << (WIFEXITED(status) ? WEXITSTATUS(status) : -1) << "\n";
        return 1;
    }
    return 0;
}

// ── CLI entry point ───────────────────────────────────────────────────────────
static void query_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " query [options]\n"
        "\n"
        "Options:\n"
        "  --db <path>          Input GBZ database file\n"
        "  --out <path>         Output GFA file\n"
        "  --sample <name>      Reference sample name (e.g. CHM13)\n"
        "  --contig <name>      Contig/chromosome name (e.g. chr1)\n"
        "  --interval <s>..<e>  0-based half-open genomic interval\n"
        "  --query-bin <path>   Path to gbz-base query binary (auto-detected if omitted)\n"
        "\n"
        "Example:\n"
        "  pva query --db graph.db --out region.gfa \\\n"
        "            --sample CHM13 --contig chr1 --interval 1000000..1001000\n";
}

int query_main(int argc, char* argv[]) {
    std::string sample, contig, interval, db, out_gfa, query_bin;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "pva query: " << a << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if      (a == "--db")        db        = next();
        else if (a == "--out")       out_gfa   = next();
        else if (a == "--sample")    sample    = next();
        else if (a == "--contig")    contig    = next();
        else if (a == "--interval")  interval  = next();
        else if (a == "--query-bin") query_bin = next();
        else if (a == "--help" || a == "-h") { query_usage(argv[0]); return 0; }
        else {
            std::cerr << "pva query: unknown option " << a << "\n";
            query_usage(argv[0]);
            return 1;
        }
    }

    if (sample.empty() || contig.empty() || interval.empty()
        || db.empty() || out_gfa.empty()) {
        query_usage(argv[0]);
        return 1;
    }

    // Parse interval string "start..end"
    auto dot = interval.find("..");
    if (dot == std::string::npos) {
        std::cerr << "pva query: --interval must be in the form start..end\n";
        return 1;
    }
    uint64_t interval_start, interval_end;
    try {
        interval_start = std::stoull(interval.substr(0, dot));
        interval_end   = std::stoull(interval.substr(dot + 2));
    } catch (...) {
        std::cerr << "pva query: invalid interval: " << interval << "\n";
        return 1;
    }
    if (interval_start >= interval_end) {
        std::cerr << "pva query: interval start must be < end\n";
        return 1;
    }

    return pva::query(db, sample, contig,
                      interval_start, interval_end,
                      out_gfa, query_bin);
}

} // namespace pva
