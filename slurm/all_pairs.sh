#!/bin/bash
# Slurm script to run centrolign as a pairwise aligner on all pairs.
# Adapted from: github.com/miramastoras/centrolign_analysis/.../all_pairs.sh
#
# Required --export variables:
#   COMBINATIONS_FILE  — path to combinations.txt (path1 \t path2, one pair per line)
#   OUTDIR             — root output directory (must contain pairwise_cigar/ and work/)
#   CENTROLIGN_BIN     — path to centrolign binary (optional; defaults to repo bin/centrolign)
#
# Submit from $OUTDIR:
#   sbatch --job-name=vntr_all_pairs \
#          --array=[1-NPAIRS]%32 \
#          --export=COMBINATIONS_FILE=...,OUTDIR=...,CENTROLIGN_BIN=... \
#          /path/to/slurm/all_pairs.sh

#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

WORKDIR="$OUTDIR/work"
CIGAR_DIR="$OUTDIR/pairwise_cigar"

mkdir -p "$CIGAR_DIR" "$WORKDIR"

# Extract paths for this array task
FASTA1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1)
FASTA2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2)

# Sample name = filename without extension (full haplotype identifier)
SAMPLE1=$(basename "$FASTA1" .fa)
SAMPLE2=$(basename "$FASTA2" .fa)

echo "sample 1: $SAMPLE1"
echo "sample 2: $SAMPLE2"
echo "fasta 1:  $FASTA1"
echo "fasta 2:  $FASTA2"
echo "out:      $CIGAR_DIR"

# Skip if already done (allows resuming failed arrays)
OUT="$CIGAR_DIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt"
if [ -f "$OUT" ]; then
    echo "Already done: $OUT"
    exit 0
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CENTROLIGN="${CENTROLIGN_BIN:-$ROOT/bin/centrolign}"

TEMP_FASTA="$WORKDIR/${SAMPLE1}_${SAMPLE2}.fa"
cat "$FASTA1" "$FASTA2" > "$TEMP_FASTA"

time "$CENTROLIGN" -v 3 --skip-calibration "$TEMP_FASTA" > "$OUT"

rm -f "$TEMP_FASTA"
