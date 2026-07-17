#!/usr/bin/env python3
"""Concatenate per-tree distance CSVs (as produced by scripts/compare_trees.R) into a
single long table with a `tree` column, for the new (Gaussian-window / Phyla / NeuralNJ)
methods. Kept separate from the existing aggregation on purpose.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Aggregate new tree-vs-ref distance CSVs.")
    ap.add_argument("--inputs", nargs="+", required=True, help="Per-tree CSVs (columns: metric,value).")
    ap.add_argument("--out", required=True, help="Output long CSV (columns: tree,metric,value).")
    args = ap.parse_args()

    frames = []
    for path in args.inputs:
        df = pd.read_csv(path)
        df.insert(0, "tree", Path(path).stem)
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {len(frames)} trees x {out['metric'].nunique()} metrics -> {args.out}")


if __name__ == "__main__":
    main()
