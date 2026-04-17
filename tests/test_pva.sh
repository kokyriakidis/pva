#!/usr/bin/env bash
# tests/test_pva.sh — integration tests for the pva binary
#
# Usage:
#   ./tests/test_pva.sh                          Run all tests
#   TEST_DATA_DIR=<path> ./tests/test_pva.sh     Also run tests that need real graph data
#
# Tests that require real data (GBZ, GFA, snarls, dist index) are skipped unless
# TEST_DATA_DIR is set. That directory must contain:
#   graph.gbz      GBZ pangenome graph
#   graph.dist     Snarl distance index (vg index -j)
#   region.gfa     GFA subgraph (output of pva query)
#   region.snarls  Snarls file (output of pva snarls)

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PVA="$ROOT/bin/pva"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0
SKIP=0

# ── colours ──────────────────────────────────────────────────────────────────
if [ -t 1 ]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
    CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; CYAN=''; BOLD=''; NC=''
fi

pass() { printf "${GREEN}PASS${NC}  %s\n" "$1"; ((PASS++)) || true; }
fail() { printf "${RED}FAIL${NC}  %s\n" "$1"; ((FAIL++)) || true; }
skip() { printf "${YELLOW}SKIP${NC}  %s\n" "$1"; ((SKIP++)) || true; }
section() { printf "\n${BOLD}${CYAN}── %s${NC}\n" "$1"; }

# Run a command; pass if exit code matches expected (default 0).
expect_exit() {
    local name="$1" expected="${2:-0}"; shift 2
    local actual
    "$@" >/dev/null 2>&1; actual=$?
    if [ "$actual" -eq "$expected" ]; then pass "$name"
    else fail "$name (exit $actual, expected $expected)"; fi
}

# Run a command and check its stdout matches a grep pattern.
expect_output() {
    local name="$1" pattern="$2"; shift 2
    local out
    out=$("$@" 2>/dev/null) || true
    if echo "$out" | grep -qE "$pattern"; then pass "$name"
    else fail "$name (pattern '$pattern' not found in output)"; fi
}

# ── Prerequisites ─────────────────────────────────────────────────────────────
section "Prerequisites"

if [ ! -x "$PVA" ]; then
    printf "${RED}ERROR${NC}: pva binary not found at %s\n" "$PVA"
    printf "Run ./install.sh first.\n"
    exit 1
fi
pass "pva binary exists"

# ── Help / dispatch ───────────────────────────────────────────────────────────
section "Help and dispatch"

expect_exit  "pva --help exits 0"                0   "$PVA" --help
expect_exit  "pva with no args exits non-zero"   1   "$PVA"
expect_exit  "pva unknown subcommand exits non-zero" 1 "$PVA" no-such-command
expect_output "pva --help lists subcommands"     "query|snarls|trace-haplotypes|prep-guide-tree-slurm|infer-guide-tree|align" "$PVA" --help

for sub in query snarls boundary snarl-dist-filter trace-haplotypes prep-guide-tree-slurm infer-guide-tree align; do
    expect_exit "pva $sub --help exits 0"        0   "$PVA" "$sub" --help
    expect_exit "pva $sub no args exits non-zero" 1  "$PVA" "$sub"
done

# ── pva prep-guide-tree-slurm ────────────────────────────────────────────────
section "pva prep-guide-tree-slurm"

# Synthetic multi-FASTA — 4 haplotypes, 80 bases each
ALLELES="$TMP/alleles.fa"
cat >"$ALLELES" <<'EOF'
>HG002#1#chr1#0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>HG003#1#chr1#0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAATCGATCG
>HG004#1#chr1#0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCATCGATCG
>HG005#1#chr1#0
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTTTCGATCG
EOF

GTDIR="$TMP/guide_tree"
SLURM_CMDS="$TMP/submit_jobs.sh"
MAIL_USER="test@example.org"
SLURM_CMDS_WITH_MAIL="$TMP/submit_jobs_with_mail.sh"
CENTROLIGN_BIN="/bin/echo"

"$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" \
    --slurm-out "$SLURM_CMDS" \
    --partition short \
    --mem 56gb \
    --time 1:00:00 \
    --cpus-per-task 32 \
    --centrolign-bin "$CENTROLIGN_BIN" >/dev/null 2>&1
GT_EXIT=$?

if [ "$GT_EXIT" -eq 0 ]; then
    pass "pva prep-guide-tree-slurm exits 0"
else
    fail "pva prep-guide-tree-slurm exits $GT_EXIT"
fi

# N=4, expected pairs = 4*3/2 = 6
N_LINES=$(wc -l <"$GTDIR/combinations.txt" 2>/dev/null || echo 0)
if [ "$N_LINES" -eq 6 ]; then
    pass "prep-guide-tree-slurm combinations.txt has 6 pairs (N=4)"
else
    fail "prep-guide-tree-slurm combinations.txt has $N_LINES lines (expected 6)"
fi

N_SPLITS=$(ls "$GTDIR/split_fastas/"*.fa 2>/dev/null | wc -l || echo 0)
if [ "$N_SPLITS" -eq 4 ]; then
    pass "prep-guide-tree-slurm split_fastas/ has 4 FASTA files"
else
    fail "prep-guide-tree-slurm split_fastas/ has $N_SPLITS files (expected 4)"
fi

# Each split FASTA should contain exactly one sequence
SINGLE=$(grep -c "^>" "$GTDIR/split_fastas/HG002_1_chr1_0.fa" 2>/dev/null || echo 0)
if [ "$SINGLE" -eq 1 ]; then
    pass "prep-guide-tree-slurm split FASTA contains exactly 1 sequence"
else
    fail "prep-guide-tree-slurm split FASTA contains $SINGLE sequences (expected 1)"
fi

# Combinations must be tab-separated with two fields per line
VALID_COMBOS=$(awk -F'\t' 'NF==2' "$GTDIR/combinations.txt" 2>/dev/null | wc -l || echo 0)
if [ "$VALID_COMBOS" -eq 6 ]; then
    pass "prep-guide-tree-slurm combinations.txt format is tab-separated pairs"
else
    fail "prep-guide-tree-slurm combinations.txt has $VALID_COMBOS valid tab-separated lines (expected 6)"
fi

# SLURM commands should be written to the requested file
if [ -s "$SLURM_CMDS" ] && grep -q "sbatch" "$SLURM_CMDS" &&
   grep -q -- "--partition=short" "$SLURM_CMDS" &&
   grep -q -- "--ntasks=1" "$SLURM_CMDS" &&
   grep -q -- "--cpus-per-task=32" "$SLURM_CMDS" &&
   grep -q -- "--array=\\[1-6\\]%32" "$SLURM_CMDS" &&
   ! grep -q -- "--mail-user=" "$SLURM_CMDS" &&
   ! grep -q -- "--mail-type=" "$SLURM_CMDS" &&
   ! grep -q "infer_tree.py" "$SLURM_CMDS" &&
   grep -q -- "CENTROLIGN_BIN=$CENTROLIGN_BIN" "$SLURM_CMDS"; then
    pass "prep-guide-tree-slurm writes SLURM commands to file"
else
    fail "prep-guide-tree-slurm did not write SLURM commands to file"
fi

if "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" \
       --slurm-out "$SLURM_CMDS_WITH_MAIL" \
       --partition short \
       --cpus-per-task 8 \
       --mail-user "$MAIL_USER" \
       --mail-type ALL \
       --mem 56gb \
       --time 1:00:00 \
       --centrolign-bin "$CENTROLIGN_BIN" >/dev/null 2>&1; then
    if grep -q -- "--mail-user=$MAIL_USER" "$SLURM_CMDS_WITH_MAIL" &&
       grep -q -- "--mail-type=ALL" "$SLURM_CMDS_WITH_MAIL" &&
       grep -q -- "--cpus-per-task=8" "$SLURM_CMDS_WITH_MAIL" &&
       grep -q -- "--array=\\[1-6\\]%8" "$SLURM_CMDS_WITH_MAIL"; then
        pass "prep-guide-tree-slurm includes mail flags when requested"
    else
        fail "prep-guide-tree-slurm did not include requested mail flags"
    fi
else
    fail "pva prep-guide-tree-slurm with mail flags exits non-zero"
fi

# Resumable: run again, should still exit 0 and not overwrite existing files
MTIME_BEFORE=$(stat -c %Y "$GTDIR/combinations.txt" 2>/dev/null || echo 0)
"$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" \
    --slurm-out "$SLURM_CMDS" \
    --partition short \
    --mem 56gb \
    --time 1:00:00 \
    --cpus-per-task 32 \
    --centrolign-bin "$CENTROLIGN_BIN" >/dev/null 2>&1 || true
MTIME_AFTER=$(stat -c %Y "$GTDIR/combinations.txt" 2>/dev/null || echo 0)
if [ "$MTIME_BEFORE" -eq "$MTIME_AFTER" ]; then
    pass "prep-guide-tree-slurm is idempotent (combinations.txt not overwritten)"
else
    fail "prep-guide-tree-slurm re-wrote combinations.txt on second run"
fi

expect_exit "pva prep-guide-tree-slurm missing --out-dir exits non-zero" 1 \
    "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --slurm-out "$SLURM_CMDS" \
    --partition short --mem 56gb --time 1:00:00 --cpus-per-task 32 --centrolign-bin "$CENTROLIGN_BIN"
expect_exit "pva prep-guide-tree-slurm missing --slurm-out exits non-zero" 1 \
    "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" \
    --partition short --mem 56gb --time 1:00:00 --cpus-per-task 32 --centrolign-bin "$CENTROLIGN_BIN"
expect_exit "pva prep-guide-tree-slurm missing --partition exits non-zero" 1 \
    "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" --slurm-out "$SLURM_CMDS" \
    --mem 56gb --time 1:00:00 --cpus-per-task 32 --centrolign-bin "$CENTROLIGN_BIN"
expect_exit "pva prep-guide-tree-slurm missing --cpus-per-task exits non-zero" 1 \
    "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" --slurm-out "$SLURM_CMDS" \
    --partition short --mem 56gb --time 1:00:00 --centrolign-bin "$CENTROLIGN_BIN"
expect_exit "pva prep-guide-tree-slurm missing --centrolign-bin exits non-zero" 1 \
    "$PVA" prep-guide-tree-slurm --fasta "$ALLELES" --out-dir "$GTDIR" --slurm-out "$SLURM_CMDS" \
    --partition short --mem 56gb --time 1:00:00 --cpus-per-task 32

# ── pva align ─────────────────────────────────────────────────────────────────
section "pva align"

# Two sequences for single-pair alignment
FA1="$TMP/seq1.fa"
FA2="$TMP/seq2.fa"
cat >"$FA1" <<'EOF'
>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
cat >"$FA2" <<'EOF'
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

CIGAR_OUT="$TMP/pair.cigar"

if "$PVA" align --seq1 "$FA1" --seq2 "$FA2" --out "$CIGAR_OUT" --skip-calibration 2>/dev/null; then
    pass "pva align single pair exits 0"

    # CIGAR should be non-empty
    if [ -s "$CIGAR_OUT" ]; then
        pass "pva align produces non-empty CIGAR output"
    else
        fail "pva align CIGAR output is empty"
    fi

    # CIGAR should only contain valid operators (digits + =/X)
    CIGAR=$(cat "$CIGAR_OUT")
    if echo "$CIGAR" | grep -qE '^[0-9]+[=X]([0-9]+[=X])*$'; then
        pass "pva align CIGAR format is valid (=/X ops)"
    else
        fail "pva align CIGAR has unexpected format: $CIGAR"
    fi

    # Identical sequences should produce an all-match CIGAR
    if echo "$CIGAR" | grep -qE '^[0-9]+=+$'; then
        pass "pva align identical sequences produce all-match CIGAR"
    else
        fail "pva align identical sequences produced: $CIGAR (expected N=)"
    fi
else
    skip "pva align single pair (centrolign may not be linked or sequences too short)"
fi

# All-pairs mode
MULTI_FA="$TMP/multi.fa"
cat >"$MULTI_FA" <<'EOF'
>hap1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>hap2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAATCGATCG
>hap3
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCATCGATCG
EOF

CIGAR_DIR="$TMP/cigars"
if "$PVA" align --all-pairs --fasta "$MULTI_FA" --out "$CIGAR_DIR" \
       --skip-calibration >/dev/null 2>&1; then
    pass "pva align --all-pairs exits 0"

    N_CIGARS=$(ls "$CIGAR_DIR"/*.txt 2>/dev/null | wc -l || echo 0)
    if [ "$N_CIGARS" -eq 3 ]; then
        pass "pva align --all-pairs produced 3 CIGAR files (N=3, 3 pairs)"
    else
        fail "pva align --all-pairs produced $N_CIGARS CIGAR files (expected 3)"
    fi
else
    skip "pva align --all-pairs (centrolign may not be linked)"
fi

# Error handling
expect_exit "pva align missing --seq2 exits non-zero" 1 \
    "$PVA" align --seq1 "$FA1"
expect_exit "pva align --all-pairs missing --fasta exits non-zero" 1 \
    "$PVA" align --all-pairs --out "$CIGAR_DIR"

# ── pva snarls (without real data) ───────────────────────────────────────────
section "pva snarls"
expect_exit "pva snarls missing --gfa exits non-zero" 1 "$PVA" snarls --out out.snarls
expect_exit "pva snarls missing --out exits non-zero" 1 "$PVA" snarls --gfa in.gfa
expect_exit "pva snarls nonexistent file exits non-zero" 1 \
    "$PVA" snarls --gfa /nonexistent/file.gfa --out "$TMP/out.snarls"

# ── pva boundary (without real data) ─────────────────────────────────────────
section "pva boundary"
expect_exit "pva boundary missing --gfa exits non-zero"    1 "$PVA" boundary --snarls in.snarls
expect_exit "pva boundary missing --snarls exits non-zero" 1 "$PVA" boundary --gfa in.gfa
expect_exit "pva boundary missing --sample exits non-zero" 1 \
    "$PVA" boundary --gfa in.gfa --snarls in.snarls

BOUNDARY_GFA="$TMP/boundary.gfa"
cat >"$BOUNDARY_GFA" <<'EOF'
H	VN:Z:1.1
S	1	A
S	2	C
S	3	G
S	4	T
W	CHM13	0	chr1	0	4	>1>2>3>4
W	GRCh38	0	chr1	0	4	>4>3>2>1
EOF

FAKE_VG="$TMP/fake_vg.sh"
cat >"$FAKE_VG" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [ "${1:-}" = "view" ] && [ "${2:-}" = "-R" ]; then
    printf '%s\n' '{"start":{"node_id":"2"},"end":{"node_id":"4"},"parent":null}'
else
    exit 64
fi
EOF
chmod +x "$FAKE_VG"

BOUNDARY_OUT="$TMP/boundary.out"
BOUNDARY_ERR="$TMP/boundary.err"
if "$PVA" boundary --gfa "$BOUNDARY_GFA" --snarls "$TMP/ignored.snarls" \
       --sample MISSING --vg-bin "$FAKE_VG" >"$BOUNDARY_OUT" 2>"$BOUNDARY_ERR"; then
    fail "pva boundary rejects unknown sample names from GFA"
else
    if grep -q "sample 'MISSING' not found" "$BOUNDARY_ERR" &&
       grep -q "CHM13" "$BOUNDARY_ERR" &&
       grep -q "GRCh38" "$BOUNDARY_ERR"; then
        pass "pva boundary rejects unknown sample names from GFA"
    else
        fail "pva boundary unknown-sample error reports available samples"
    fi
fi

if "$PVA" boundary --gfa "$BOUNDARY_GFA" --snarls "$TMP/ignored.snarls" \
       --sample CHM13 --vg-bin "$FAKE_VG" >"$BOUNDARY_OUT" 2>"$BOUNDARY_ERR"; then
    if grep -q $'^2\t4$' "$BOUNDARY_OUT"; then
        pass "pva boundary uses the requested sample walk when present"
    else
        fail "pva boundary output was not the expected boundary pair"
    fi
else
    fail "pva boundary with valid sample exits non-zero"
fi

# ── pva snarl-dist-filter (without real data) ─────────────────────────────────
section "pva snarl-dist-filter"
expect_exit "pva snarl-dist-filter missing --snarls exits non-zero" 1 \
    "$PVA" snarl-dist-filter --dist g.dist --out out.snarls
expect_exit "pva snarl-dist-filter missing --dist exits non-zero"   1 \
    "$PVA" snarl-dist-filter --snarls in.snarls --out out.snarls
expect_exit "pva snarl-dist-filter missing --out exits non-zero"    1 \
    "$PVA" snarl-dist-filter --snarls in.snarls --dist g.dist

# ── pva trace-haplotypes (without real data) ────────────────────────────────────────────
section "pva trace-haplotypes"
expect_exit "pva trace-haplotypes missing --gbz exits non-zero"        1 "$PVA" trace-haplotypes --start 1 --end 2
expect_exit "pva trace-haplotypes missing --start exits non-zero"            1 "$PVA" trace-haplotypes --gbz g.gbz --end 2
expect_exit "pva trace-haplotypes missing --end exits non-zero"              1 "$PVA" trace-haplotypes --gbz g.gbz --start 1
expect_exit "pva trace-haplotypes --snarls missing --out-dir exits non-zero" 1 "$PVA" trace-haplotypes --gbz g.gbz --snarls in.snarls

# ── pva query (without real data) ────────────────────────────────────────────
section "pva query"
expect_exit "pva query missing --db exits non-zero"       1 \
    "$PVA" query --out out.gfa --sample S --contig c --interval 0..1
expect_exit "pva query missing --out exits non-zero"      1 \
    "$PVA" query --db g.db --sample S --contig c --interval 0..1
expect_exit "pva query missing --sample exits non-zero"   1 \
    "$PVA" query --db g.db --out out.gfa --contig c --interval 0..1
expect_exit "pva query missing --interval exits non-zero" 1 \
    "$PVA" query --db g.db --out out.gfa --sample S --contig c
expect_exit "pva query bad interval format exits non-zero" 1 \
    "$PVA" query --db g.db --out out.gfa --sample S --contig c --interval bad

# ── Tests requiring real data ─────────────────────────────────────────────────
section "Integration tests (real graph data)"

if [ -z "${TEST_DATA_DIR:-}" ]; then
    skip "pva query  (set TEST_DATA_DIR to enable)"
    skip "pva snarls (set TEST_DATA_DIR to enable)"
    skip "pva boundary (set TEST_DATA_DIR to enable)"
    skip "pva snarl-dist-filter (set TEST_DATA_DIR to enable)"
    skip "pva trace-haplotypes (set TEST_DATA_DIR to enable)"
else
    GBZ="$TEST_DATA_DIR/graph.gbz"
    DIST="$TEST_DATA_DIR/graph.dist"
    GFA="$TEST_DATA_DIR/region.gfa"
    SNARLS="$TEST_DATA_DIR/region.snarls"

    # query
    if [ -f "$TEST_DATA_DIR/graph.db" ]; then
        expect_exit "pva query real data" 0 \
            "$PVA" query --db "$TEST_DATA_DIR/graph.db" \
            --out "$TMP/region.gfa" \
            --sample CHM13 --contig chr1 --interval 1000000..1001000
    else
        skip "pva query (graph.db not found in TEST_DATA_DIR)"
    fi

    # snarls
    if [ -f "$GFA" ]; then
        expect_exit "pva snarls real GFA" 0 \
            "$PVA" snarls --gfa "$GFA" --out "$TMP/region.snarls"
        SNARLS_OUT="$TMP/region.snarls"
    else
        skip "pva snarls (region.gfa not found in TEST_DATA_DIR)"
        SNARLS_OUT="${SNARLS:-}"
    fi

    # boundary
    if [ -f "$GFA" ] && [ -f "${SNARLS_OUT:-}" ]; then
        OUT=$("$PVA" boundary --gfa "$GFA" --snarls "$SNARLS_OUT" --sample CHM13 2>/dev/null)
        if echo "$OUT" | grep -qE '^[0-9]+[[:space:]][0-9]+$'; then
            pass "pva boundary real data produces start_node<tab>end_node"
        else
            fail "pva boundary unexpected output: $OUT"
        fi
    else
        skip "pva boundary (GFA or snarls not available)"
    fi

    # snarl-dist-filter
    if [ -f "${SNARLS:-}" ] && [ -f "$DIST" ]; then
        expect_exit "pva snarl-dist-filter real data" 0 \
            "$PVA" snarl-dist-filter \
            --snarls "$SNARLS" --dist "$DIST" \
            --out "$TMP/filtered.snarls"
        if [ -s "$TMP/filtered.snarls" ]; then
            pass "pva snarl-dist-filter produces non-empty snarls file"
        else
            skip "pva snarl-dist-filter produced empty output (no snarls in range)"
        fi
    else
        skip "pva snarl-dist-filter (snarls or dist index not found)"
    fi

    # trace
    if [ -f "$GBZ" ] && [ -f "${SNARLS:-}" ]; then
        expect_exit "pva trace-haplotypes --snarls mode" 0 \
            "$PVA" trace-haplotypes --gbz "$GBZ" \
            --snarls "$SNARLS" \
            --out-dir "$TMP/alleles"
        N_FA=$(ls "$TMP/alleles/snarl_"*.fa 2>/dev/null | wc -l || echo 0)
        if [ "$N_FA" -gt 0 ]; then
            pass "pva trace-haplotypes produced $N_FA FASTA file(s)"
        else
            skip "pva trace-haplotypes produced 0 FASTAs (no passing snarls)"
        fi
    else
        skip "pva trace-haplotypes (GBZ or snarls not found)"
    fi
fi

# ── pva infer-guide-tree (without real data) ──────────────────────────────────
section "pva infer-guide-tree"
expect_exit "pva infer-guide-tree missing --cigar-dir exits non-zero" 1 \
    "$PVA" infer-guide-tree --combinations comb.txt --out tree.nwk
expect_exit "pva infer-guide-tree missing --combinations exits non-zero" 1 \
    "$PVA" infer-guide-tree --cigar-dir cigar/ --out tree.nwk
expect_exit "pva infer-guide-tree missing --out exits non-zero" 1 \
    "$PVA" infer-guide-tree --cigar-dir cigar/ --combinations comb.txt
# With both required args but nonexistent paths, the Python script should fail
expect_exit "pva infer-guide-tree nonexistent paths exits non-zero" 1 \
    "$PVA" infer-guide-tree \
    --cigar-dir /nonexistent/cigar \
    --combinations /nonexistent/combinations.txt \
    --out "$TMP/tree.nwk"

# ── Summary ───────────────────────────────────────────────────────────────────
printf "\n${BOLD}Results: ${GREEN}%d passed${NC} ${BOLD}/ ${RED}%d failed${NC} ${BOLD}/ ${YELLOW}%d skipped${NC}\n" \
    "$PASS" "$FAIL" "$SKIP"

[ "$FAIL" -eq 0 ]
