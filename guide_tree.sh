#!/usr/bin/env bash
# guide_tree.sh
# Usage: guide_tree.sh <haplotypes.fa> <outdir>
#
# Prepares all-vs-all centrolign run, following the same pipeline as:
# github.com/miramastoras/centrolign_analysis/.../submit_all_pairs.md
#
# Steps:
#   1. Split multi-FASTA into per-haplotype FASTA files
#   2. Build fasta_list.txt (one absolute path per line)
#   3. Generate all-vs-all combinations.txt (path1 \t path2)
#   4. Print sbatch commands (split into chunks if needed)
#
# After submitting all SLURM jobs, run:
#   python3 scripts/infer_tree.py <outdir>/pairwise_cigar <outdir>/combinations.txt \
#       > <outdir>/guide_tree.nwk

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <haplotypes.fa> <outdir>"
    exit 1
fi

FASTA=$(realpath "$1")
OUTDIR=$(realpath "$2")

ROOT="$(cd "$(dirname "$0")" && pwd)"
SLURM_SCRIPT="$ROOT/slurm/all_pairs.sh"

SPLIT_DIR="$OUTDIR/split_fastas"
FASTA_LIST="$OUTDIR/fasta_list.txt"
COMBINATIONS="$OUTDIR/combinations.txt"
CIGAR_DIR="$OUTDIR/pairwise_cigar"
LOG_DIR="$OUTDIR/logs"

mkdir -p "$SPLIT_DIR" "$CIGAR_DIR" "$LOG_DIR"

# ── 1. Split multi-FASTA into per-haplotype files ────────────────────────────
echo "Splitting haplotypes..."

python3 - "$FASTA" "$SPLIT_DIR" <<'PYEOF'
import sys, os, re

fasta_file = sys.argv[1]
out_dir    = sys.argv[2]

def safe(label):
    return re.sub(r'[^0-9a-zA-Z.\-]', '_', label)

label = seq_parts = None
count = 0
with open(fasta_file) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if label is not None and seq_parts:
                path = os.path.join(out_dir, safe(label) + '.fa')
                with open(path, 'w') as out:
                    out.write(f'>{label}\n{"".join(seq_parts)}\n')
                count += 1
            label = line[1:]
            seq_parts = []
        elif label is not None:
            seq_parts.append(line)
    if label is not None and seq_parts:
        path = os.path.join(out_dir, safe(label) + '.fa')
        with open(path, 'w') as out:
            out.write(f'>{label}\n{"".join(seq_parts)}\n')
        count += 1

print(f"{count} haplotypes written to {out_dir}")
PYEOF

# ── 2. Build fasta_list.txt ───────────────────────────────────────────────────
echo "Building fasta list..."
> "$FASTA_LIST"
for f in "$SPLIT_DIR"/*.fa; do
    realpath "$f"
done | sort >> "$FASTA_LIST"

N=$(wc -l < "$FASTA_LIST")
echo "Haplotypes: $N"

# ── 3. Generate all-vs-all combinations (lexicographic, no duplicate pairs) ──
echo "Generating all-vs-all combinations..."

while read -r s1; do
    while read -r s2; do
        [[ "$s1" < "$s2" ]] && printf '%s\t%s\n' "$s1" "$s2"
    done < "$FASTA_LIST"
done < "$FASTA_LIST" > "$COMBINATIONS"

NPAIRS=$(wc -l < "$COMBINATIONS")
echo "Total pairs: $NPAIRS"

# ── Sanity check: no empty FASTAs ────────────────────────────────────────────
EMPTY=$(find "$SPLIT_DIR" -type f -empty | wc -l)
if [ "$EMPTY" -gt 0 ]; then
    echo "WARNING: $EMPTY empty FASTA file(s) in $SPLIT_DIR"
    find "$SPLIT_DIR" -type f -empty
fi

# ── 4. Print sbatch commands (chunked to respect array job limit) ─────────────
MAX_ARRAY=30000   # SLURM array limit on most clusters
CONCURRENT=128    # max simultaneously running tasks per batch

echo ""
echo "=== Setup complete ==="
echo ""
echo "Submit the following SLURM array job(s):"
echo ""
echo "  cd $OUTDIR"
echo ""

start=1
batch=1
while [ "$start" -le "$NPAIRS" ]; do
    end=$(( start + MAX_ARRAY - 1 ))
    [ "$end" -gt "$NPAIRS" ] && end="$NPAIRS"
    echo "  # Batch $batch (tasks ${start}-${end})"
    echo "  sbatch \\"
    echo "      --job-name=vntr_all_pairs_${batch} \\"
    echo "      --array=[${start}-${end}]%${CONCURRENT} \\"
    echo "      --export=COMBINATIONS_FILE=${COMBINATIONS},OUTDIR=${OUTDIR} \\"
    echo "      $SLURM_SCRIPT"
    echo ""
    start=$(( end + 1 ))
    (( batch++ ))
done

echo "After all jobs complete, infer the guide tree:"
echo "  python3 $ROOT/scripts/infer_tree.py $CIGAR_DIR $COMBINATIONS > $OUTDIR/guide_tree.nwk"
