#!/usr/bin/env python3
"""
infer_tree.py
Usage: infer_tree.py <cigar_dir> <combinations.txt> [> guide_tree.nwk]

Reads pairwise centrolign CIGAR files, builds a distance matrix:

    distance = 1 - (2 * matches) / (ref_len + query_len)

Infers a neighbor-joining tree (scikit-bio) and prints Newick to stdout.

combinations.txt: tab-separated path1  path2  (one pair per line)
CIGAR files:      <cigar_dir>/pairwise_cigar_<sample1>_<sample2>.txt
"""

import sys
import os
import re

try:
    import skbio
except ImportError:
    sys.exit("ERROR: scikit-bio required.  Install: pip install scikit-bio")


def parse_cigar(text):
    return [(m.group(2), int(m.group(1)))
            for m in re.finditer(r'(\d+)([HSMIDX=])', text)]


def cigar_to_dist(cigar):
    query_len = ref_len = matches = 0
    for op, length in cigar:
        if op in ('M', 'X', '='):
            query_len += length
            ref_len   += length
            if op != 'X':
                matches += length
        elif op == 'D':
            ref_len   += length
        elif op in ('I', 'H', 'S'):
            query_len += length
    denom = ref_len + query_len
    return 1.0 if denom == 0 else 1.0 - (2.0 * matches) / denom


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: infer_tree.py <cigar_dir> <combinations.txt>")

    cigar_dir  = sys.argv[1]
    combos     = sys.argv[2]

    pairs = []
    samples_seen = set()
    with open(combos) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            # paths in combinations.txt; sample name = basename without .fa
            name1 = os.path.basename(parts[0]).removesuffix('.fa')
            name2 = os.path.basename(parts[1]).removesuffix('.fa')
            pairs.append((name1, name2))
            samples_seen.update([name1, name2])

    samples = sorted(samples_seen)
    n       = len(samples)
    idx     = {s: i for i, s in enumerate(samples)}

    D = [[0.0] * n for _ in range(n)]

    missing = []
    for name1, name2 in pairs:
        cigar_file = os.path.join(cigar_dir, f"pairwise_cigar_{name1}_{name2}.txt")
        if not os.path.isfile(cigar_file):
            # Try reversed order (shouldn't happen since we use lexicographic pairs)
            cigar_file = os.path.join(cigar_dir, f"pairwise_cigar_{name2}_{name1}.txt")
        if not os.path.isfile(cigar_file):
            missing.append((name1, name2))
            continue

        with open(cigar_file) as f:
            text = f.read().strip()

        dist = 1.0 if not text else cigar_to_dist(parse_cigar(text))
        i, j = idx[name1], idx[name2]
        D[i][j] = D[j][i] = dist

    if missing:
        print(f"WARNING: {len(missing)} CIGAR file(s) missing — "
              "treating as max distance (1.0)", file=sys.stderr)
        for name1, name2 in missing:
            i, j = idx[name1], idx[name2]
            D[i][j] = D[j][i] = 1.0

    dist_mat = skbio.DistanceMatrix(D, samples)
    tree     = skbio.tree.nj(dist_mat)
    tree     = tree.root_at_midpoint()
    print(tree)


if __name__ == '__main__':
    main()
