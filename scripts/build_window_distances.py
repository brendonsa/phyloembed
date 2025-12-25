#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from compute_distance_matrix_phylip import write_phylip_relaxed


def compute_windows(L: int, window_size: int, overlap: int):
    if window_size <= 0:
        raise ValueError("window_size must be > 0")
    if overlap < 0:
        raise ValueError("overlap must be >= 0")
    if overlap >= window_size:
        raise ValueError("overlap must be < window_size")

    step = window_size - overlap
    windows = []
    start = 0

    while start + window_size <= L:
        windows.append((start, start + window_size))
        start += step

    if start < L:
        windows.append((start, L))  # tail window

    return windows


def main():
    ap = argparse.ArgumentParser(
        description="Write per-window PHYLIP distance files and an aggregated joined.phy from concat embeddings NPZ."
    )
    ap.add_argument("--embeddings", required=True, help="Input .npz (keys: labels, embeddings[N,L,D])")
    ap.add_argument("--outdir", required=True, help="Output directory (writes distance_XXXX.phy + joined.phy + windows.json)")
    ap.add_argument("--window-size", type=int, required=True, help="Window size in residues")
    ap.add_argument("--overlap", type=int, default=0, help="Overlap in residues (default: 0)")
    ap.add_argument("--metric", type=str, default="cosine", help="pdist metric (default: cosine)")
    ap.add_argument("--digits", type=int, default=4, help="Zero pad for window index (default: 4)")
    ap.add_argument("--reduce", choices=["mean", "sum"], default="mean", help="Aggregate over windows (default: mean)")
    ap.add_argument("--prefix", default="distance_", help="Per-window filename prefix (default: distance_)")
    ap.add_argument(
        "--joined-only",
        action="store_true",
        help="Only write joined.phy (skip per-window .phy files)",
    )
    args = ap.parse_args()

    data = np.load(args.embeddings, allow_pickle=True)
    if "labels" not in data or "embeddings" not in data:
        raise ValueError("NPZ must contain keys: 'labels' and 'embeddings'.")

    labels = data["labels"]          # (N,)
    emb = data["embeddings"]         # (N,L,D)
    if emb.ndim != 3:
        raise ValueError(f"Expected embeddings shape (N,L,D), got {emb.shape}")

    N, L, D = emb.shape
    windows = compute_windows(L, args.window_size, args.overlap)
    if not windows:
        raise ValueError("No windows produced (check window_size/overlap).")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    acc = None
    files = []

    for i, (s, e) in enumerate(windows):
        W = emb[:, s:e, :].mean(axis=1)                  # (N, D)
        dist = squareform(pdist(W, metric=args.metric))  # (N, N)

        if not args.joined_only:
            fname = f"{args.prefix}{i:0{args.digits}d}.phy"
            fpath = outdir / fname
            write_phylip_relaxed(dist, labels, str(fpath))
            files.append(fname)

        if acc is None:
            acc = dist
        else:
            acc += dist

    if args.reduce == "mean":
        acc = acc / float(len(windows))

    joined_path = outdir / "joined.phy"
    write_phylip_relaxed(acc, labels, str(joined_path))

    meta = {
        "embeddings_file": os.path.abspath(args.embeddings),
        "outdir": os.path.abspath(str(outdir)),
        "n_seq": int(N),
        "seq_len": int(L),
        "emb_dim": int(D),
        "metric": args.metric,
        "window_size": int(args.window_size),
        "overlap": int(args.overlap),
        "reduce": args.reduce,
        "n_windows": int(len(windows)),
        "windows": [{"start": int(s), "end": int(e)} for (s, e) in windows],
        "per_window_files": files,
        "joined_file": "joined.phy",
        "joined_only": bool(args.joined_only),
    }
    with open(outdir / "windows.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote {len(windows)} window PHYLIP files + joined.phy to {outdir}")


if __name__ == "__main__":
    main()
