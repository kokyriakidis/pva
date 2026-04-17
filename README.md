# pangenome-vntr-alignments

A toolkit for extracting VNTR-spanning haplotypes from pangenome graphs and
preparing them for downstream comparison. `pva` can query a genomic interval
from a GBZ database, find snarls and boundary nodes in the resulting subgraph,
trace haplotype sequences across a single snarl or an entire region, and then
prepare all-vs-all centrolign alignments for guide-tree inference.

> Development status: `pva` is still under active development and should not
> be used for production workflows or relied on for stable results yet.

## Installation

```bash
./install.sh
```

Installs/builds the required dependencies, builds `bin/pva` (+ `bin/libpva.so`),
and prepares the bundled tests.

All required external tools (`vg`, `gbz-base`, `centrolign`) are downloaded or
built into `bin/` and auto-detected at runtime — no PATH setup needed. The
`--vg-bin` / `--query-bin` override flags exist only for non-standard
installations where you want to point to a different binary.

---

## Commands

### `pva query` — Extract a GFA subgraph for a genomic interval

**Required inputs:**

| Flag | Description |
|---|---|
| `--db <path>` | GBZ database file (`.db`) produced by `gbz2db` |
| `--out <path>` | Output GFA file |
| `--sample <name>` | Reference sample name (e.g. `CHM13`) |
| `--contig <name>` | Contig / chromosome name (e.g. `chr1`) |
| `--interval <s>..<e>` | 0-based half-open genomic interval |

```bash
pva query \
  --db graph.db \
  --out region.gfa \
  --sample CHM13 \
  --contig chr1 \
  --interval 1000000..1001000
```

---

### `pva snarls` — Find snarls in a GFA subgraph

**Required inputs:**

| Flag | Description |
|---|---|
| `--gfa <path>` | Input GFA file (output of `pva query`) |
| `--out <path>` | Output snarls file |

```bash
pva snarls \
  --gfa region.gfa \
  --out region.snarls
```

---

### `pva snarl-dist-filter` — Filter snarls by minimum distance between boundary nodes

Queries `vg distance` for every top-level snarl and keeps those whose minimum
distance between boundary nodes falls within the given range. Requires the
distance index built by `vg index -j`.

**Required inputs:**

| Flag | Description |
|---|---|
| `--snarls <path>` | Snarls file (output of `pva snarls`) |
| `--dist <path>` | Snarl distance index (`.dist`) built by `vg index -j` |
| `--out <path>` | Output TSV file |

**Optional:**

| Flag | Description |
|---|---|
| `--min-dist <N>` | Keep snarls with min distance ≥ N (default: `0`) |
| `--max-dist <N>` | Keep snarls with min distance ≤ N (default: unlimited) |

**Output:** a filtered binary snarls file, directly consumable by `pva trace-haplotypes --snarls`.

```bash
pva snarl-dist-filter \
  --snarls region.snarls \
  --dist graph.dist \
  --out filtered.snarls \
  --min-dist 100 \
  --max-dist 50000
```

---

### `pva boundary` — Find outermost boundary nodes of top-level snarls

Scans the reference sample's walk (W line) in the GFA to determine which
boundary nodes are genomically first and last across all top-level snarls
combined. Useful when you want to trace all haplotypes across the whole region
as one unit rather than per-snarl.

**Required inputs:**

| Flag | Description |
|---|---|
| `--gfa <path>` | Input GFA file |
| `--snarls <path>` | Snarls file (output of `pva snarls`) |
| `--sample <name>` | Reference sample for walk ordering; must exist in a GFA `W` line |

**Output:** `start_node\tend_node` written to stdout.

If the requested sample is not present in the GFA, `pva boundary` exits with an
error and reports the available sample names found in `W` lines.

```bash
pva boundary \
  --gfa region.gfa \
  --snarls region.snarls \
  --sample CHM13
```

---

### `pva trace-haplotypes` — Trace haplotype sequences between two nodes

Extracts all haplotype sequences that span from a start node to an end node
in the GBZ graph, following all GBWT-indexed paths via DFS. Works across any
range — a single snarl, a multi-snarl region, or an arbitrary pair of nodes.

The r-index (~15–50× faster locate()) is built automatically on first run and
cached as `<graph.gbz>.ri`. On subsequent runs it is loaded from disk. Use
`--ri` to specify a custom path if the cache lives elsewhere.

**Optional flag available in all modes:**

| Flag | Description |
|---|---|
| `--ri <path>` | Pre-built r-index file (default: `<graph.gbz>.ri`, built and cached on first run) |

#### Use case 1 — Per-snarl tracing (normal pipeline)

Iterates over all snarls in a filtered snarls file, producing one multi-FASTA
per snarl. Each file contains all haplotypes that span that snarl's boundary
nodes. This is the output fed to `pva prep-guide-tree-slurm` and `pva align`.

| Flag | Description |
|---|---|
| `--gbz <path>` | Input GBZ graph file |
| `--snarls <path>` | Filtered snarls file (output of `pva snarl-dist-filter`) |
| `--out-dir <path>` | Output directory; one file per snarl named `snarl_<start>_<end>.fa` |

```bash
pva trace-haplotypes --gbz graph.gbz --snarls filtered.snarls --out-dir snarls/
# On first run, builds and caches graph.gbz.ri automatically.
# On subsequent runs, loads it from disk.
# To use a pre-built r-index at a custom path:
pva trace-haplotypes --gbz graph.gbz --ri /data/graph.gbz.ri --snarls filtered.snarls --out-dir snarls/
```

Output `snarls/` directory:
```
snarls/snarl_12345_12360.fa   ← all haplotypes spanning snarl 12345→12360
snarls/snarl_23000_23100.fa   ← all haplotypes spanning snarl 23000→23100
...
```

Each file contains one sequence per haplotype:
```
>HG002#1#chr1#0
ATCGATCG...
>HG003#1#chr1#0
ATCGATCG...
```

#### Use case 2 — Whole-region tracing (spans multiple snarls)

Traces all haplotypes between any two arbitrary nodes — the DFS crosses snarl
boundaries freely. Use `pva boundary` to find the outermost nodes of the whole
region automatically.

| Flag | Description |
|---|---|
| `--gbz <path>` | Input GBZ graph file |
| `--start <node>` | Start boundary node ID |
| `--end <node>` | End boundary node ID |
| `--out <path>` | Output multi-FASTA (default: stdout) |

```bash
# Find the outermost boundary nodes of the whole region
BOUNDS=$(pva boundary --gfa region.gfa --snarls region.snarls --sample CHM13)
START=$(echo "$BOUNDS" | cut -f1)
END=$(echo "$BOUNDS"   | cut -f2)

# Trace all haplotypes spanning the entire region (across all snarls)
pva trace-haplotypes --gbz graph.gbz --start "$START" --end "$END" --out whole_region.fa

# With a pre-built r-index:
pva trace-haplotypes --gbz graph.gbz --ri /data/graph.gbz.ri --start "$START" --end "$END" --out whole_region.fa
```

For multiple custom ranges, loop `pva trace-haplotypes` in shell:

```bash
while IFS=$'\t' read -r start end; do
    pva trace-haplotypes --gbz graph.gbz --start "$start" --end "$end" \
              --out "region_${start}_${end}.fa"
done < ranges.tsv
```

---

### `pva prep-guide-tree-slurm` — Prepare all-vs-all alignment for guide tree inference

Splits a multi-FASTA into per-haplotype files, generates all N×(N−1)/2
pairwise combinations, and writes the SLURM `sbatch` commands to a file.
After jobs complete, run `pva infer-guide-tree` to produce a Newick guide tree.

This step is intended for cluster execution. The number of pairwise alignments
grows quadratically with the number of haplotypes, so even when individual
centrolign runs are fast, whole-region guide tree preparation can still require
substantial total compute time.

**Required inputs:**

| Flag | Description |
|---|---|
| `--fasta <path>` | Multi-FASTA of haplotype sequences (output of `pva trace-haplotypes`) |
| `--out-dir <path>` | Output directory for split FASTAs, `combinations.txt`, logs, and pairwise CIGARs |
| `--slurm-out <path>` | Output file that will contain the generated SLURM submission commands |
| `--partition <name>` | SLURM partition/queue name |
| `--mem <spec>` | SLURM memory request, e.g. `56gb` |
| `--time <spec>` | SLURM time limit, e.g. `1:00:00` |
| `--centrolign-bin <path>` | Path to the `centrolign` binary that worker jobs should run |
| `--cpus-per-task <N>` | SLURM CPUs per task and array parallelism cap |

**Optional:**

| Flag | Description |
|---|---|
| `--nodes <N>` | SLURM node count (default: `1`) |
| `--log-pattern <path>` | SLURM log path pattern relative to `--out-dir` (default: `logs/array_job_%A_task_%a.log`) |
| `--job-name-prefix <name>` | Prefix for generated SLURM job names (default: `vntr_all_pairs`) |
| `--mail-user <addr>` | Email address for SLURM notifications |
| `--mail-type <type>` | SLURM mail type, e.g. `ALL`, `END`, or `FAIL` |
| `--max-array <N>` | SLURM array job size limit (default: `30000`) |

**Creates inside `--out-dir`:**

```text
guide_tree/
├── combinations.txt     # all pairwise FASTA combinations
├── fasta_list.txt       # absolute paths to split FASTAs
├── logs/                # SLURM stdout/stderr logs
├── pairwise_cigar/      # pairwise_cigar_*.txt outputs from SLURM jobs
└── split_fastas/        # one FASTA per haplotype
```

`--slurm-out` is written separately and contains the generated `sbatch`
commands. Run `pva infer-guide-tree` later as a separate step after the SLURM
jobs finish (see below).

```bash
pva prep-guide-tree-slurm \
  --fasta alleles.fa \
  --out-dir guide_tree/ \
  --slurm-out guide_tree/submit_jobs.sh \
  --partition short \
  --mem 56gb \
  --time 1:00:00 \
  --cpus-per-task 32 \
  --centrolign-bin /path/to/centrolign
```

---

### `pva infer-guide-tree` — Infer NJ guide tree from pairwise CIGAR files

Runs `scripts/infer_tree.py` to build a distance matrix from the completed
pairwise CIGAR files and infer a neighbor-joining tree using scikit-bio.

**Required inputs:**

| Flag | Description |
|---|---|
| `--cigar-dir <path>` | Directory containing `pairwise_cigar_*.txt` files produced by the completed SLURM jobs, typically `<out-dir>/pairwise_cigar` |
| `--combinations <path>` | `combinations.txt` created by `pva prep-guide-tree-slurm` |
| `--out <path>` | Output Newick file |

**Optional:**

| Flag | Description |
|---|---|
| `--script <path>` | Path to `infer_tree.py` (auto-detected if omitted) |

```bash
pva infer-guide-tree \
  --cigar-dir guide_tree/pairwise_cigar \
  --combinations guide_tree/combinations.txt \
  --out guide_tree/guide_tree.nwk
```

#### Guide-tree workflows

#### Use case 1 — Build guide tree from interval coordinates

```bash
# 1. Extract a regional subgraph from genomic coordinates
pva query \
  --db graph.db \
  --out region.gfa \
  --sample CHM13 \
  --contig chr1 \
  --interval 1000000..2000000

# 2. Find snarls in that region
pva snarls --gfa region.gfa --out region.snarls

# 3. Find the outermost boundary nodes of the whole region
BOUNDS=$(pva boundary --gfa region.gfa --snarls region.snarls --sample CHM13)
START=$(echo "$BOUNDS" | cut -f1)
END=$(echo "$BOUNDS" | cut -f2)

# 4. Trace all haplotypes spanning the full region
pva trace-haplotypes \
  --gbz graph.gbz \
  --start "$START" \
  --end "$END" \
  --out whole_region.fa

# 5. Prepare SLURM array jobs for guide tree construction
pva prep-guide-tree-slurm \
  --fasta whole_region.fa \
  --out-dir guide_tree/ \
  --slurm-out guide_tree/submit_jobs.sh \
  --partition short \
  --mem 56gb \
  --time 1:00:00 \
  --cpus-per-task 32 \
  --centrolign-bin /path/to/centrolign

# 6. After the SLURM jobs finish, infer the Newick tree
pva infer-guide-tree \
  --cigar-dir guide_tree/pairwise_cigar \
  --combinations guide_tree/combinations.txt \
  --out guide_tree/guide_tree.nwk
```

#### Use case 2 — Build guide tree from known start and end nodes

```bash
# 1. Trace all haplotypes spanning the chosen node interval
pva trace-haplotypes \
  --gbz graph.gbz \
  --start 12345 \
  --end 67890 \
  --out region.fa

# 2. Prepare SLURM array jobs for guide tree construction
pva prep-guide-tree-slurm \
  --fasta region.fa \
  --out-dir guide_tree/ \
  --slurm-out guide_tree/submit_jobs.sh \
  --partition short \
  --mem 56gb \
  --time 1:00:00 \
  --cpus-per-task 32 \
  --centrolign-bin /path/to/centrolign

# 3. After the SLURM jobs finish, infer the Newick tree
pva infer-guide-tree \
  --cigar-dir guide_tree/pairwise_cigar \
  --combinations guide_tree/combinations.txt \
  --out guide_tree/guide_tree.nwk
```

---

### `pva align` — Pairwise or all-vs-all alignment via centrolign

Runs centrolign alignment in-memory (no subprocess). Outputs an explicit CIGAR
string using `=`/`X` operators.

**Single pair:**

| Flag | Description |
|---|---|
| `--seq1 <path>` | First input FASTA |
| `--seq2 <path>` | Second input FASTA |

**All-vs-all:**

| Flag | Description |
|---|---|
| `--all-pairs` | Run all N×(N−1)/2 pairs from a multi-FASTA |
| `--fasta <path>` | Multi-FASTA input |
| `--out <path>` | Output directory for CIGAR files |

**Optional (both modes):**

| Flag | Description |
|---|---|
| `--out <path>` | Output CIGAR file (single-pair, default: stdout) |
| `--skip-calibration` | Skip anchor score calibration (faster; recommended for guide tree) |
| `-v <level>` | Verbosity: `0`=silent `1`=minimal `2`=basic `3`=verbose (default: `0`) |

```bash
# Single pair
pva align --seq1 HG002.fa --seq2 HG003.fa

# All-vs-all
pva align --all-pairs --fasta alleles.fa --out cigar_dir/ --skip-calibration
```

---

## Binary overrides

After `./install.sh` all tools are bundled in `bin/` and auto-detected — these
flags are not needed for normal use. They exist only if you want to substitute
a specific version of an external tool.

| Flag | Applies to | Default location |
|---|---|---|
| `--vg-bin <path>` | `snarls`, `boundary`, `snarl-dist-filter` | `bin/vg` |
| `--query-bin <path>` | `query` | `bin/gbz-base/target/release/query` |

---

## File types reference

| Extension | Description |
|---|---|
| `.db` | GBZ database (from `gbz2db`) |
| `.gbz` | Pangenome graph (GBWT index + GBWTGraph) |
| `.gfa` | Graph Fragment Assembly subgraph |
| `.snarls` | Binary snarls file (from `vg snarls`) |
| `.dist` | Snarl distance index (from `vg index -j`) |
| `.ri` | GBWT r-index cache (auto-created by `pva trace-haplotypes` on first run; specify with `--ri`) |
| `.fa` | FASTA sequence file |
| `.nwk` | Newick guide tree |
