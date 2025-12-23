#!/usr/bin/env python3
import argparse
from Bio import Phylo
from Bio.Phylo.Consensus import get_support, majority_consensus


def strip_internal_labels(tree):
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None  # drop inner labels that would clash with support


def main(base_tree_path, tree_paths, out_supported, cutoff):
    if not tree_paths:
        raise SystemExit("No bootstrap trees provided via --trees.")
    base = Phylo.read(base_tree_path, "newick")
    strip_internal_labels(base)

    
    trees = []
    for p in tree_paths:
        print(p)
        t = Phylo.read(p, "newick")
        strip_internal_labels(t)
        trees.append(t)

    supported = get_support(base, trees)
    Phylo.write(supported, out_supported, "newick")

    # consensus = majority_consensus(trees, cutoff=cutoff)
    # Phylo.write(consensus, out_consensus, "newick")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Annotate base tree with supports and write consensus (no globbing).")
    ap.add_argument("--base-tree", required=True, help="Path to base tree (.nwk)")
    ap.add_argument("--trees", nargs="+", required=True, help="List of bootstrap tree files (.nwk)")
    ap.add_argument("--out-supported", required=True, help="Output path for base-with-supports (.nwk)")
    # ap.add_argument("--out-consensus", required=True, help="Output path for consensus (.nwk)")
    ap.add_argument("--cutoff", type=float, default=0.5, help="Majority-rule cutoff (default 0.5)")
    args = ap.parse_args()
    main(args.base_tree, args.trees, args.out_supported, args.cutoff)
