#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from compute_distance_matrix_phylip import write_phylip_relaxed


def main():
    ap = argparse.ArgumentParser(
        description="Sample weighted site subsets (bucket-filling) and write PHYLIP distance matrices per subset + joined.phy."
    )
    ap.add_argument("--embeddings", required=True, help="Input .npz with keys: labels, embeddings (N,L,D)")
    ap.add_argument("--weights", required=True, help="Input .npy site weights of length L")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--n-subsets", type=int, required=True, help="Number of subset distance matrices to produce")
    ap.add_argument("--coverage", type=float, default=0.5, help="Coverage fraction of total positive weight (default: 0.5)")
    ap.add_argument("--metric", type=str, default="cosine", help="pdist metric (default: cosine)")
    ap.add_argument("--include-zero-weights", action="store_true", help="Allow selecting w<=0 sites (still adds 0 to coverage)")
    ap.add_argument("--random-state", type=int, default=1337, help="RNG seed (default: 1337)")
    ap.add_argument("--prefix", default="subset_", help="Per-subset filename prefix (default: subset_)")
    ap.add_argument("--digits", type=int, default=4, help="Zero-pad digits for subset index (default: 4)")
    ap.add_argument("--reduce", choices=["mean", "sum"], default="mean", help="How to aggregate subsets into joined.phy (default: mean)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = np.load(args.embeddings, allow_pickle=True)
    if "labels" not in data or "embeddings" not in data:
        raise ValueError("NPZ must contain keys: 'labels' and 'embeddings'.")

    labels = data["labels"]          # (N,)
    emb = data["embeddings"]         # (N,L,D)
    if emb.ndim != 3:
        raise ValueError(f"Expected embeddings shape (N,L,D). Got {emb.shape}")

    N, L, D = emb.shape

    site_weights = np.asarray(np.load(args.weights), dtype=float)
    if site_weights.ndim != 1:
        raise ValueError("weights must be a 1D array")
    if len(site_weights) != L:
        raise ValueError(f"weights length {len(site_weights)} != embeddings L {L}")

    pos_mask = site_weights > 0.0
    total_pos_weight = float(site_weights[pos_mask].sum())
    if total_pos_weight <= 0:
        raise ValueError("No positive site weights found.")

    if args.include_zero_weights:
        candidate_indices = np.arange(L, dtype=int)
    else:
        candidate_indices = np.where(pos_mask)[0].astype(int)

    if len(candidate_indices) == 0:
        raise ValueError("No candidate sites to sample (check include_zero_weights/weights).")

    target_weight = float(args.coverage) * total_pos_weight
    rng = np.random.default_rng(args.random_state)

    subset_files = []
    subsets = []  # list of dicts: selected sites, cum weight, count
    joined_acc = None
    produced = 0

    for sidx in range(int(args.n_subsets)):
        # shuffled permutation of candidate sites; each site used at most once in this subset
        indices = candidate_indices.copy()
        rng.shuffle(indices)

        selected = []
        cum_w = 0.0

        for idx in indices:
            w = float(site_weights[idx])
            if (w <= 0.0) and (not args.include_zero_weights):
                continue
            selected.append(int(idx))
            cum_w += max(w, 0.0)
            if cum_w >= target_weight:
                break

        if not selected:
            continue

        # aggregate: mean of per-site distance matrices (unweighted)
        agg = None
        for idx in selected:
            W = emb[:, idx, :]                          # (N, D) site embeddings
            dist = squareform(pdist(W, metric=args.metric))  # (N, N)
            if agg is None:
                agg = dist
            else:
                agg += dist
        agg = agg / float(len(selected))

        fname = f"{args.prefix}{sidx:0{args.digits}d}.phy"
        write_phylip_relaxed(agg, labels, str(outdir / fname))
        subset_files.append(fname)

        if joined_acc is None:
            joined_acc = agg.copy()
        else:
            joined_acc += agg

        subsets.append(
            {
                "subset_index": int(sidx),
                "n_selected_sites": int(len(selected)),
                "cum_positive_weight": float(cum_w),
                "target_weight": float(target_weight),
                "selected_sites": selected,
                "file": fname,
            }
        )
        produced += 1

    if produced == 0:
        raise ValueError("Produced 0 subsets (check coverage/include_zero_weights/weights).")

    if args.reduce == "mean":
        joined_acc = joined_acc / float(produced)

    joined_path = outdir / "joined.phy"
    write_phylip_relaxed(joined_acc, labels, str(joined_path))

    meta = {
        "embeddings_file": os.path.abspath(args.embeddings),
        "weights_file": os.path.abspath(args.weights),
        "outdir": os.path.abspath(str(outdir)),
        "n_seq": int(N),
        "seq_len": int(L),
        "emb_dim": int(D),
        "metric": args.metric,
        "coverage": float(args.coverage),
        "include_zero_weights": bool(args.include_zero_weights),
        "random_state": int(args.random_state),
        "requested_subsets": int(args.n_subsets),
        "produced_subsets": int(produced),
        "reduce": args.reduce,
        "subset_files": subset_files,
        "subsets": subsets,
        "joined_file": "joined.phy",
    }
    with open(outdir / "subsets.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote {produced} subset PHYLIP files + joined.phy to {outdir}")


if __name__ == "__main__":
    main()
