

# Convert a GBZ graph graph.gbz into a database graph.db

# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/gbz-base/target/release/gbz2db \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.gbz \
# --output /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.db
#--overwrite


# gbz-base repository at: https://github.com/jltsiren/gbz-base
# vg repository at: https://github.com/vgteam/vg

# # Generate the snarls file from the graph
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg snarls \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.gbz \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.snarls


# # Snarl-based query.
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/gbz-base/target/release/query \
# --context 0 \
# --between 9380487:9380493 \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/kolf2.1j-dg-sample.db \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/kolf2.1j-dg-sample-9380487-9380493.gfa




#!/usr/bin/env bash
set -euo pipefail

GBZ=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.gbz
DB=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.db
GFA=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.gfa
SNARLS=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.snarls
OUT=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region_alleles.fa

VG=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg
QUERY=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/gbz-base/target/release/query
TRACE=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/trace_haplotypes

# 1. Query by genomic interval → subgraph GFA
"$QUERY" --sample CHM13 --contig chr1 --interval 1000000..1001000 "$DB" > "$GFA"

# 2. Find snarls in the subgraph
"$VG" snarls "$GFA" > "$SNARLS"

# 3. Find genomic start/end boundary nodes.
# Uses the CHM13 reference walk from the GFA to determine order —
# correct even when node IDs are not assigned in genomic order.
read START_NODE END_NODE < <(
  "$VG" view -R "$SNARLS" | \
  jq -r 'select(.parent == null) | .start.node_id, .end.node_id' | \
  awk -v gfa="$GFA" '
    { boundary[$1] = 1 }
    END {
      while ((getline line < gfa) > 0) {
        if (line ~ /^W\tCHM13\t/) {
          n = split(line, f, "\t")
          walk = f[7]
          while (match(walk, /[><][0-9]+/)) {
            node = substr(walk, RSTART+1, RLENGTH-1)
            walk = substr(walk, RSTART+RLENGTH)
            if (node in boundary) {
              if (start == "") start = node
              end = node
            }
          }
          break
        }
      }
      print start, end
    }
  '
)

echo "Genomic start boundary node: $START_NODE"
echo "Genomic end boundary node:   $END_NODE"

if [ -z "$START_NODE" ] || [ -z "$END_NODE" ] || [ "$START_NODE" = "$END_NODE" ]; then
  echo "ERROR: could not determine valid boundary nodes"
  exit 1
fi

# 4. Trace haplotypes through the snarl and emit labeled FASTA
"$TRACE" "$GBZ" "$START_NODE" "$END_NODE" > "$OUT"

SEQ_COUNT=$(grep -c "^>" "$OUT" 2>/dev/null || echo 0)
echo "Sequences extracted: $SEQ_COUNT"






# # 1. Query by genomic interval → subgraph
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/gbz-base/target/release/query \
# --sample CHM13 \
# --contig chr1 \
# --interval 1000000..1001000 \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.db \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.gfa

# # 2. Find snarls in that subgraph
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg snarls \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.gfa \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.snarls

# # 3. Find leftmost and rightmost boundary nodes across all top-level snarls
# # node_id is a string in vg JSON output; tonumber converts for numeric sort.
# # Top-level snarls have no "parent" field (comes through as null in jq).
# GFA=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.gfa
# SNARLS=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.snarls

# # Collect all top-level snarl boundary node IDs into a set.
# BOUNDARY_NODES=$(
#   /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg view -R "$SNARLS" | \
#   jq -r 'select(.parent == null) | .start.node_id, .end.node_id'
# )

# # Find the CHM13 reference walk in the GFA and extract nodes in genomic order.
# # W lines: sample hap contig start end walk — we take the first CHM13 walk.
# REF_WALK=$(grep -m1 "^W	CHM13	" "$GFA" | awk '{print $7}')

# # Walk through the reference walk in order, emit boundary nodes as encountered.
# # The first boundary node seen = genomic start; the last = genomic end.
# START_NODE=""
# END_NODE=""
# while IFS= read -r node; do
#     if echo "$BOUNDARY_NODES" | grep -qx "$node"; then
#         [ -z "$START_NODE" ] && START_NODE="$node"
#         END_NODE="$node"
#     fi
# done < <(echo "$REF_WALK" | grep -oP '(?<=[><])\d+')

# echo "Genomic start boundary node: $START_NODE"
# echo "Genomic end boundary node:   $END_NODE"

# # Verify both nodes appear as S lines in the GFA.
# if ! grep -qP "^S\t${START_NODE}\t" "$GFA"; then
#   echo "ERROR: start node $START_NODE not in GFA — boundary selection may be wrong"
#   exit 1
# fi
# if ! grep -qP "^S\t${END_NODE}\t" "$GFA"; then
#   echo "ERROR: end node $END_NODE not in GFA — boundary selection may be wrong"
#   exit 1
# fi

# # Count paths visiting start node (sanity check).
# SPANNING=$(grep -cP "(>|<)${START_NODE}(>|<)" "$GFA" 2>/dev/null || echo 0)
# echo "Paths visiting start node in GFA: $SPANNING"
# if [ "$SPANNING" -eq 0 ]; then
#   echo "WARNING: no paths visit start node $START_NODE — output will be empty"
# fi

# # 4. Trace haplotypes through the snarl and emit labeled FASTA
# OUT=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region_alleles.fa
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/trace_haplotypes \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/hprc-v2.1-mc-chm13-eval.gbz \
# "$START_NODE" "$END_NODE" \
# > "$OUT"

# # Report how many sequences were produced.
# SEQ_COUNT=$(grep -c "^>" "$OUT" 2>/dev/null || echo 0)
# echo "Sequences extracted: $SEQ_COUNT"
# if [ "$SEQ_COUNT" -eq 0 ]; then
#   echo "WARNING: no output — check snarl pairs above and consider using a specific pair"
# fi





# # 3. Extract top-level snarl boundary nodes
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg view -R /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.snarls | \
#   jq 'select(.parent == null or .parent.start.node_id == 0) | 
#     {start: .start.node_id, end: .end.node_id}'

# # Build filter_paths (only recompiles if source is newer than binary)
# FILTER_SRC=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/filter_paths.cpp
# FILTER_BIN=/home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/filter_paths
# if [ ! -f "$FILTER_BIN" ] || [ "$FILTER_SRC" -nt "$FILTER_BIN" ]; then
#     g++ -O2 -std=c++17 -o "$FILTER_BIN" "$FILTER_SRC"
# fi

# # 4. Filter paths that don't span both boundary nodes
# START_NODE=12345
# END_NODE=12346
# "$FILTER_BIN" \
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.gfa \
# "$START_NODE" "$END_NODE" \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.filtered.gfa

# # 5. Extract path sequences
# /home/kokyriakidis/Downloads/pangenome-vntr-alignments/bin/vg paths \
# --extract-fasta \
# -x /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/test_region.filtered.gfa \
# > /home/kokyriakidis/Downloads/pangenome-vntr-alignments/data/region_sequences.fa
