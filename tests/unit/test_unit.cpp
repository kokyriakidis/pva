// tests/unit/test_unit.cpp — unit tests for pva utility functions
//
// Tests pva_utils.h in isolation (no heavy deps: no gbwt, no centrolign, no vg).
// Compiled standalone by install.sh into tests/unit/test_unit.
//
// Run: ./tests/unit/test_unit

#include "../../pva_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <sys/stat.h>

// ── Minimal test framework ────────────────────────────────────────────────────

static int s_pass = 0, s_fail = 0;

static void test(bool cond, const char* name) {
    if (cond) {
        std::cout << "\033[0;32mPASS\033[0m  " << name << "\n";
        ++s_pass;
    } else {
        std::cout << "\033[0;31mFAIL\033[0m  " << name << "\n";
        ++s_fail;
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

static std::string tmp_path(const std::string& suffix) {
    const char* base = std::getenv("TMPDIR");
    if (!base) base = "/tmp";
    return std::string(base) + "/pva_test_" + suffix;
}

static void write_file(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    return std::string(std::istreambuf_iterator<char>(f),
                       std::istreambuf_iterator<char>());
}

// ── Tests: find_bin ───────────────────────────────────────────────────────────

static void test_find_bin() {
    std::cout << "\n── find_bin ─────────────────────────────────────────────\n";

    // hint takes priority
    test(pva::utils::find_bin("/my/custom/tool", "relative", "fallback")
         == "/my/custom/tool",
         "find_bin: hint is returned as-is");

    // empty hint falls through to PATH fallback
    std::string r = pva::utils::find_bin("", "nonexistent_rel_binary_xyz", "sh");
    test(r == "sh",
         "find_bin: falls back to PATH name when relative path not found");

    // /proc/self/exe is readable — self_dir detection works
    // (relative path exists would return it; here we just check it doesn't crash)
    std::string r2 = pva::utils::find_bin("", "nonexistent_binary_abc", "cat");
    test(!r2.empty(), "find_bin: always returns a non-empty string");
}

// ── Tests: run_capture ────────────────────────────────────────────────────────

static void test_run_capture() {
    std::cout << "\n── run_capture ──────────────────────────────────────────\n";

    // echo hello
    std::string out;
    int rc = pva::utils::run_capture({"echo", "hello"}, out);
    test(rc == 0, "run_capture: echo exits 0");
    test(out == "hello\n", "run_capture: echo output is 'hello\\n'");

    // capture multiple lines
    std::string out2;
    pva::utils::run_capture({"printf", "a\\nb\\nc\\n"}, out2);
    test(out2 == "a\nb\nc\n", "run_capture: multi-line output");

    // non-zero exit is propagated
    std::string out3;
    int rc3 = pva::utils::run_capture({"sh", "-c", "exit 42"}, out3);
    test(rc3 == 42, "run_capture: non-zero exit code propagated");

    // stdout captured, not stderr
    std::string out4;
    pva::utils::run_capture({"sh", "-c", "echo stdout; echo stderr >&2"}, out4);
    test(out4 == "stdout\n", "run_capture: captures stdout only");
}

// ── Tests: run_to_file ────────────────────────────────────────────────────────

static void test_run_to_file() {
    std::cout << "\n── run_to_file ──────────────────────────────────────────\n";

    std::string path = tmp_path("run_to_file.txt");

    int rc = pva::utils::run_to_file({"echo", "written"}, path);
    test(rc == 0, "run_to_file: exits 0");
    test(read_file(path) == "written\n", "run_to_file: file content correct");

    // Overwrites existing file
    pva::utils::run_to_file({"echo", "overwritten"}, path);
    test(read_file(path) == "overwritten\n", "run_to_file: overwrites existing file");

    // Non-zero exit
    int rc2 = pva::utils::run_to_file({"sh", "-c", "exit 7"}, path);
    test(rc2 == 7, "run_to_file: non-zero exit code propagated");

    std::remove(path.c_str());
}

// ── Tests: run_pipe_to_file ───────────────────────────────────────────────────

static void test_run_pipe_to_file() {
    std::cout << "\n── run_pipe_to_file ─────────────────────────────────────\n";

    std::string path = tmp_path("pipe_to_file.txt");

    // cat echoes stdin to stdout
    std::string input = "line1\nline2\nline3\n";
    int rc = pva::utils::run_pipe_to_file({"cat"}, input, path);
    test(rc == 0, "run_pipe_to_file: cat exits 0");
    test(read_file(path) == input, "run_pipe_to_file: output matches input");

    // Empty input
    pva::utils::run_pipe_to_file({"cat"}, "", path);
    test(read_file(path).empty(), "run_pipe_to_file: empty input produces empty file");

    // Large input (> 65 kB to exercise the write loop)
    std::string large(200000, 'A');
    large.back() = '\n';
    pva::utils::run_pipe_to_file({"cat"}, large, path);
    test(read_file(path) == large, "run_pipe_to_file: handles large input (200 kB)");

    // Transform via tr
    std::string upper_input = "hello world\n";
    pva::utils::run_pipe_to_file({"tr", "a-z", "A-Z"}, upper_input, path);
    test(read_file(path) == "HELLO WORLD\n",
         "run_pipe_to_file: transform via tr");

    std::remove(path.c_str());
}

// ── Tests: snarl JSON parsing (via regex-based logic) ─────────────────────────
// We reproduce the same parsing logic used in pva_snarl_dist_filter.cpp
// so we can test it in isolation without the full binary.

static bool line_is_top_level(const std::string& line) {
    auto pos = line.find("\"parent\"");
    if (pos == std::string::npos) return true;
    auto colon = line.find(':', pos);
    if (colon == std::string::npos) return true;
    size_t v = colon + 1;
    while (v < line.size() && line[v] == ' ') ++v;
    return line.substr(v, 4) == "null";
}

static void test_snarl_parsing() {
    std::cout << "\n── snarl JSON top-level detection ──────────────────────\n";

    // No parent field → top-level
    test(line_is_top_level(R"({"start":{"node_id":"1"},"end":{"node_id":"2"}})"),
         "snarl with no parent field is top-level");

    // parent: null → top-level
    test(line_is_top_level(R"({"start":{"node_id":"3"},"end":{"node_id":"4"},"parent":null})"),
         "snarl with parent:null is top-level");

    // parent: {} → NOT top-level
    test(!line_is_top_level(R"({"start":{"node_id":"5"},"end":{"node_id":"6"},"parent":{"node_id":"1"}})"),
         "snarl with parent object is NOT top-level");

    // parent with extra whitespace
    test(line_is_top_level(R"({"parent":  null})"),
         "snarl with parent:  null (extra space) is top-level");
}

// ── Tests: FASTA format (smoke test via pva_guide_tree logic) ─────────────────
// We test the safe_name function behaviour expected by guide-tree.
// The function replaces non-alphanumeric chars (except . and -) with _.

static std::string safe_name_ref(const std::string& label) {
    std::string r;
    for (char c : label) {
        bool ok = std::isalnum((unsigned char)c) || c == '.' || c == '-';
        r += ok ? c : '_';
    }
    return r;
}

static void test_safe_name() {
    std::cout << "\n── safe_name (guide-tree FASTA naming) ─────────────────\n";

    test(safe_name_ref("HG002#1#chr1#0") == "HG002_1_chr1_0",
         "safe_name: # replaced with _");
    test(safe_name_ref("sample.hap-1") == "sample.hap-1",
         "safe_name: . and - preserved");
    test(safe_name_ref("CHM13 v2.0") == "CHM13_v2.0",
         "safe_name: spaces and dots handled");
    test(safe_name_ref("abc123") == "abc123",
         "safe_name: alphanumeric unchanged");
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "pva unit tests\n";

    test_find_bin();
    test_run_capture();
    test_run_to_file();
    test_run_pipe_to_file();
    test_snarl_parsing();
    test_safe_name();

    std::cout << "\n";
    if (s_fail == 0)
        std::cout << "\033[0;32mAll " << s_pass << " tests passed.\033[0m\n";
    else
        std::cout << "\033[0;31m" << s_fail << " test(s) FAILED\033[0m  ("
                  << s_pass << " passed)\n";

    return s_fail > 0 ? 1 : 0;
}
