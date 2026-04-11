// filter_paths.cpp
// Usage: filter_paths <graph.gfa> <start_node> <end_node>
//
// Reads a GFA file and writes to stdout only the path/walk lines (P and W)
// that visit both boundary nodes. S, L, and H lines are always passed through.
//
// Uses mmap for zero-copy input and a large output buffer for throughput.

#include <cstdio>
#include <cstring>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

// ── output buffer ─────────────────────────────────────────────────────────────

static constexpr size_t OUT_BUF = 1 << 22; // 4 MB
static char   out_buf[OUT_BUF];
static size_t out_pos = 0;

static void flush_out() {
    if (out_pos) { fwrite(out_buf, 1, out_pos, stdout); out_pos = 0; }
}

static void write_out(const char* data, size_t len) {
    while (len) {
        size_t space = OUT_BUF - out_pos;
        size_t chunk = len < space ? len : space;
        memcpy(out_buf + out_pos, data, chunk);
        out_pos += chunk; data += chunk; len -= chunk;
        if (out_pos == OUT_BUF) flush_out();
    }
}

// ── field extraction ──────────────────────────────────────────────────────────

// Return the nth tab-delimited field as a string_view (zero-copy).
static std::string_view get_field(std::string_view line, int n) {
    size_t start = 0;
    for (int i = 0; i < n; ++i) {
        size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) return {};
        start = tab + 1;
    }
    size_t end = line.find('\t', start);
    return line.substr(start, end == std::string_view::npos ? end : end - start);
}

// ── node scanning ─────────────────────────────────────────────────────────────

// Scan a W-line walk string (e.g. ">12344>12345<12346") for two target node IDs.
// Short-circuits as soon as both are found — no allocations.
static bool walk_has_both(std::string_view walk,
                           std::string_view start, std::string_view end) {
    bool found_start = false, found_end = false;
    size_t i = 0;
    while (i < walk.size()) {
        if (walk[i] == '>' || walk[i] == '<') {
            ++i;
            size_t j = i;
            while (j < walk.size() && walk[j] != '>' && walk[j] != '<') ++j;
            std::string_view node = walk.substr(i, j - i);
            if (!found_start && node == start) found_start = true;
            if (!found_end   && node == end)   found_end   = true;
            if (found_start && found_end) return true;
            i = j;
        } else {
            ++i;
        }
    }
    return false;
}

// Scan a P-line segment field (e.g. "12344+,12345+,12346-") for two target node IDs.
// Short-circuits as soon as both are found.
static bool path_has_both(std::string_view segs,
                           std::string_view start, std::string_view end) {
    bool found_start = false, found_end = false;
    size_t i = 0;
    while (i <= segs.size()) {
        size_t j = segs.find(',', i);
        if (j == std::string_view::npos) j = segs.size();
        std::string_view seg = segs.substr(i, j - i);
        if (!seg.empty() && (seg.back() == '+' || seg.back() == '-'))
            seg.remove_suffix(1);
        if (!found_start && seg == start) found_start = true;
        if (!found_end   && seg == end)   found_end   = true;
        if (found_start && found_end) return true;
        i = j + 1;
    }
    return false;
}

// ── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <graph.gfa> <start_node> <end_node>\n", argv[0]);
        return 1;
    }

    const std::string_view start_node = argv[2];
    const std::string_view end_node   = argv[3];

    // mmap the input file — no per-line copies
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { perror(argv[1]); return 1; }

    struct stat st;
    fstat(fd, &st);
    size_t file_size = static_cast<size_t>(st.st_size);

    const char* data = static_cast<const char*>(
        mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
    close(fd);
    if (data == MAP_FAILED) { perror("mmap"); return 1; }

    // Hint to kernel: sequential access pattern
    madvise(const_cast<char*>(data), file_size, MADV_SEQUENTIAL);

    const char* p   = data;
    const char* eof = data + file_size;

    while (p < eof) {
        const char* nl       = static_cast<const char*>(memchr(p, '\n', eof - p));
        const char* line_end = nl ? nl : eof;
        std::string_view line(p, line_end - p);
        p = nl ? nl + 1 : eof;

        if (line.empty()) { write_out("\n", 1); continue; }

        char type = line[0];

        if (type == 'W') {
            std::string_view walk = get_field(line, 6);
            // Drop malformed lines (fewer than 7 fields) silently
            if (!walk.empty() && walk_has_both(walk, start_node, end_node)) {
                write_out(line.data(), line.size());
                write_out("\n", 1);
            }
        } else if (type == 'P') {
            std::string_view segs = get_field(line, 2);
            // Drop malformed lines (fewer than 3 fields) silently
            if (!segs.empty() && path_has_both(segs, start_node, end_node)) {
                write_out(line.data(), line.size());
                write_out("\n", 1);
            }
        } else {
            // S, L, H — always keep
            write_out(line.data(), line.size());
            write_out("\n", 1);
        }
    }

    flush_out();
    munmap(const_cast<char*>(data), file_size);
    return 0;
}
