#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from compute_distance_matrix_phylip import write_phylip_relaxed


def aggregate_by_coverage_sorted_from_embeddings(
    embeddings_nld: np.ndarray,
    site_weights: np.ndarray,
    coverage: float,
    metric: str,
    direction: str,
    include_zero_weights: bool,
    normalize: str,
):
    # embeddings_nld: (N, L, D)
    if embeddings_nld.ndim != 3:
        raise ValueError(f"Expected embeddings shape (N,L,D), got {embeddings_nld.shape}")
    N, L, D = embeddings_nld.shape

    site_weights = np.asarray(site_weights, dtype=float)
    if site_weights.shape != (L,):
        raise ValueError(f"weights length {site_weights.shape} != L {L}")

    pos_mask = site_weights > 0.0
    total_pos_weight = float(site_weights[pos_mask].sum())
    if total_pos_weight <= 0:
        raise ValueError("No positive site weights found.")

    target_weight = float(coverage) * total_pos_weight

    if include_zero_weights:
        candidate = np.arange(L, dtype=int)
    else:
        candidate = np.where(pos_mask)[0].astype(int)

    if direction == "high":
        order = np.argsort(site_weights[candidate])[::-1]
    elif direction == "low":
        order = np.argsort(site_weights[candidate])
    else:
        raise ValueError("direction must be 'high' or 'low'")

    selected = []
    cum_w = 0.0
    for j in order:
        idx = int(candidate[j])
        w = float(site_weights[idx])
        if (w <= 0.0) and (not include_zero_weights):
            continue
        selected.append(idx)
        cum_w += max(w, 0.0)
        if cum_w >= target_weight:
            break

    if not selected:
        raise ValueError("No sites selected (check coverage/include_zero_weights/weights).")

    agg = None
    if normalize == "weights":
        den = 0.0
        for idx in selected:
            w = max(float(site_weights[idx]), 0.0)
            W = embeddings_nld[:, idx, :]                    # (N, D)
            dist = squareform(pdist(W, metric=metric))       # (N, N)
            agg = (w * dist) if agg is None else (agg + w * dist)
            den += w
        if den == 0.0:
            raise ValueError("Normalization denominator is zero (all selected weights <= 0).")
        agg = agg / den
    elif normalize == "count":
        for idx in selected:
            W = embeddings_nld[:, idx, :]
            dist = squareform(pdist(W, metric=metric))
            agg = dist if agg is None else (agg + dist)
        agg = agg / float(len(selected))
    else:
        raise ValueError("normalize must be 'weights' or 'count'")

    meta = {
        "coverage": float(coverage),
        "target_weight": float(target_weight),
        "cum_positive_weight": float(cum_w),
        "direction": direction,
        "include_zero_weights": bool(include_zero_weights),
        "normalize": normalize,
        "metric": metric,
        "n_selected_sites": int(len(selected)),
        "selected_sites": selected,
    }
    return agg, meta


def main():
    ap = argparse.ArgumentParser(
        description="Deterministically aggregate per-site distances by sorted weight coverage; write one PHYLIP file."
    )
    ap.add_argument("--embeddings", required=True, help="NPZ with keys: labels, embeddings (N,L,D)")
    ap.add_argument("--weights", required=True, help="NPY site weights length L")
    ap.add_argument("--out", required=True, help="Output PHYLIP path")
    ap.add_argument("--meta", required=True, help="Output JSON metadata path")
    ap.add_argument("--coverage", type=float, required=True, help="0..1")
    ap.add_argument("--metric", default="cosine")
    ap.add_argument("--direction", choices=["high", "low"], default="high")
    ap.add_argument("--normalize", choices=["weights", "count"], default="weights")
    ap.add_argument("--include-zero", action="store_true")
    args = ap.parse_args()

    data = np.load(args.embeddings, allow_pickle=True)
    if "labels" not in data or "embeddings" not in data:
        raise ValueError("Embeddings NPZ must contain keys: labels, embeddings")
    labels = data["labels"]
    emb = data["embeddings"]

    w = np.load(args.weights)

    dist, meta = aggregate_by_coverage_sorted_from_embeddings(
        embeddings_nld=emb,
        site_weights=w,
        coverage=float(args.coverage),
        metric=args.metric,
        direction=args.direction,
        include_zero_weights=bool(args.include_zero),
        normalize=args.normalize,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_phylip_relaxed(dist, labels, str(out_path))

    meta_path = Path(args.meta)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_full = {
        "embeddings": os.path.abspath(args.embeddings),
        "weights": os.path.abspath(args.weights),
        "out_phylip": os.path.abspath(str(out_path)),
        **meta,
    }
    with open(meta_path, "w") as f:
        json.dump(meta_full, f, indent=2)


if __name__ == "__main__":
    main()
