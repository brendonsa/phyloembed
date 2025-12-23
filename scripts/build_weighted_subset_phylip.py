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
        description="Sample weighted site subsets and write per-subset PHYLIP distance matrices (no tree building)."
    )
    ap.add_argument("--embeddings", required=True, help="NPZ with keys: labels, embeddings (N,L,D)")
    ap.add_argument("--weights", required=True, help="NPY site weights length L")
    ap.add_argument("--outdir", required=True, help="Output directory for distance_XXXX.phy")
    ap.add_argument("--n-subsets", type=int, required=True)
    ap.add_argument("--coverage", type=float, required=True, help="0..1")
    ap.add_argument("--metric", default="cosine")
    ap.add_argument("--include-zero", action="store_true")
    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--digits", type=int, default=4)
    ap.add_argument("--prefix", default="distance_")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = np.load(args.embeddings, allow_pickle=True)
    if "labels" not in data or "embeddings" not in data:
        raise ValueError("Embeddings NPZ must contain keys: labels, embeddings")

    labels = data["labels"]          # (N,)
    emb = data["embeddings"]         # (N,L,D)
    if emb.ndim != 3:
        raise ValueError(f"Expected embeddings shape (N,L,D); got {emb.shape}")

    N, L, D = emb.shape
    site_weights = np.asarray(np.load(args.weights), dtype=float)
    if site_weights.shape != (L,):
        raise ValueError(f"weights length {site_weights.shape} != L {L}")

    pos_mask = site_weights > 0.0
    total_pos = float(site_weights[pos_mask].sum())
    if total_pos <= 0:
        raise ValueError("No positive site weights found.")
    target_weight = float(args.coverage) * total_pos

    if args.include_zero:
        candidate_indices = np.arange(L, dtype=int)
    else:
        candidate_indices = np.where(pos_mask)[0].astype(int)

    rng = np.random.default_rng(args.seed)

    produced = 0
    subsets_meta = []
    phylip_files = []

    for sidx in range(int(args.n_subsets)):
        indices = candidate_indices.copy()
        rng.shuffle(indices)

        selected = []
        cum_w = 0.0
        for idx in indices:
            w = float(site_weights[idx])
            if (w <= 0.0) and (not args.include_zero):
                continue
            selected.append(int(idx))
            cum_w += max(w, 0.0)
            if cum_w >= target_weight:
                break

        if not selected:
            continue

        # aggregate distance = mean of per-site distances (unweighted)
        agg = None
        for site in selected:
            W = emb[:, site, :]                         # (N,D)
            dist = squareform(pdist(W, metric=args.metric))  # (N,N)
            agg = dist if agg is None else (agg + dist)
        agg = agg / float(len(selected))

        fname = f"{args.prefix}{sidx:0{args.digits}d}.phy"
        write_phylip_relaxed(agg, labels, str(outdir / fname))

        phylip_files.append(fname)
        subsets_meta.append(
            {
                "subset_index": int(sidx),
                "n_selected_sites": int(len(selected)),
                "cum_positive_weight": float(cum_w),
                "target_weight": float(target_weight),
                "selected_sites": selected,
                "phylip_file": fname,
            }
        )
        produced += 1

    if produced == 0:
        raise ValueError("Produced 0 subsets (check coverage/include_zero/weights).")

    meta = {
        "embeddings": os.path.abspath(args.embeddings),
        "weights": os.path.abspath(args.weights),
        "n_seq": int(N),
        "seq_len": int(L),
        "emb_dim": int(D),
        "metric": args.metric,
        "include_zero": bool(args.include_zero),
        "seed": int(args.seed),
        "coverage": float(args.coverage),
        "n_requested": int(args.n_subsets),
        "n_produced": int(produced),
        "phylip_files": phylip_files,
        "subsets": subsets_meta,
    }
    with open(outdir / "subsets.json", "w") as f:
        json.dump(meta, f, indent=2)

    # convenience “done” marker
    (outdir / "DONE").write_text("ok\n")


if __name__ == "__main__":
    main()