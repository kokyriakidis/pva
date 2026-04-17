// pva_main.cpp — CLI entry point for the pva binary
//
// Usage: pva <subcommand> [options]
//
// Subcommands:
//   query       Extract a GFA subgraph for a genomic interval
//   snarls      Find snarls in a GFA subgraph
//   boundary    Find outermost boundary nodes of top-level snarls
//   trace       Trace haplotype sequences spanning a snarl
//   guide-tree  Prepare all-vs-all centrolign run for guide tree inference

#include "pva.h"
#include <iostream>

static void usage(std::ostream& os, const char* prog) {
    os <<
        "pva — Pangenome VNTR Aligner\n"
        "\n"
        "Usage: " << prog << " <subcommand> [options]\n"
        "\n"
        "Subcommands:\n"
        "  query       Extract a GFA subgraph for a genomic interval\n"
        "  snarls      Find snarls in a GFA subgraph\n"
        "  boundary    Find outermost boundary nodes of top-level snarls\n"
        "  trace-haplotypes  Trace haplotype sequences spanning a snarl\n"
        "  prep-guide-tree-slurm  Prepare all-vs-all centrolign run + write SLURM commands\n"
        "  infer-guide-tree        Infer NJ guide tree from pairwise CIGAR files\n"
        "  snarl-dist-filter Filter snarls by min distance between boundary nodes\n"
        "  align             Pairwise or all-vs-all alignment via centrolign\n"
        "\n"
        "Full pipeline:\n"
        "  pva query             --db graph.db --out region.gfa --sample CHM13 --contig chr1 --interval 1e6..1.001e6\n"
        "  pva snarls            --gfa region.gfa --out region.snarls\n"
        "  pva snarl-dist-filter --snarls region.snarls --dist graph.dist --min-dist 100 --max-dist 50000 > regions.tsv\n"
        "  pva trace-haplotypes  --gbz graph.gbz --snarls filtered.snarls --out-dir alleles/\n"
        "  pva prep-guide-tree-slurm --fasta alleles.fa --out-dir guide_tree/ --slurm-out guide_tree/submit_jobs.sh --partition short --mem 56gb --time 1:00:00 --cpus-per-task 32 --centrolign-bin /path/to/centrolign\n"
        "  pva infer-guide-tree        --cigar-dir guide_tree/pairwise_cigar --combinations guide_tree/combinations.txt --out guide_tree/guide_tree.nwk\n"
        "\n"
        "Run `" << prog << " <subcommand> --help` for subcommand-specific options.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(std::cerr, argv[0]); return 1; }

    std::string cmd = argv[1];

    if (cmd == "query")      return pva::query_main     (argc - 1, argv + 1);
    if (cmd == "snarls")     return pva::snarls_main    (argc - 1, argv + 1);
    if (cmd == "boundary")   return pva::boundary_main  (argc - 1, argv + 1);
    if (cmd == "trace-haplotypes") return pva::trace_main(argc - 1, argv + 1);
    if (cmd == "prep-guide-tree-slurm") return pva::prep_guide_tree_slurm_main(argc - 1, argv + 1);
    if (cmd == "infer-guide-tree")      return pva::infer_guide_tree_main  (argc - 1, argv + 1);
    if (cmd == "snarl-dist-filter") return pva::snarl_dist_filter_main(argc - 1, argv + 1);
    if (cmd == "align")       return pva::align_main        (argc - 1, argv + 1);

    if (cmd == "--help" || cmd == "-h") { usage(std::cout, argv[0]); return 0; }

    std::cerr << "pva: unknown subcommand '" << cmd << "'\n\n";
    usage(std::cerr, argv[0]);
    return 1;
}
