def column_change_counts(seqs):
    """
    seqs: list of equal-length strings (aligned sequences)
    returns: list with per-column change count = (#unique chars) - 1
    """
    # sanity check
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs), "All sequences must have same length"

    change_counts = []
    for col in zip(*seqs):              # each col is a tuple of residues
        if not col:
            continue
        k = len(set(col))
        change_counts.append(max(0, k - 1))
    return change_counts


import numpy as np

def aggregate_by_coverage_sorted(
    dist_matrices,
    site_weights,
    coverage=0.5,
    direction="high",          # "high" = most informative first, "low" = least informative first
    include_zero_weights=False,
    normalize="count",       # "weights" or "count"
):
    """
    Build ONE aggregated distance matrix by deterministically selecting sites
    until the cumulative positive weight reaches `coverage * total_positive_weight`.

    dist_matrices: list/array of (n_seq x n_seq) arrays, one per site
    site_weights:  1D array-like, length = n_sites
    coverage:      fraction of total positive weight to cover (0..1)
    direction:     "high" or "low" (sort order by site_weights)
    include_zero_weights: if False, ignore w<=0 sites entirely
    normalize:     "weights" -> weighted mean by w; "count" -> mean over selected sites
    """
    site_weights = np.asarray(site_weights, dtype=float)
    n_sites = len(dist_matrices)
    if n_sites != len(site_weights):
        raise ValueError("weights and dist_matrices length mismatch")

    pos_mask = site_weights > 0.0
    total_pos_weight = site_weights[pos_mask].sum()
    if total_pos_weight <= 0:
        raise ValueError("No positive site weights found.")

    target_weight = coverage * total_pos_weight

    if include_zero_weights:
        candidate = np.arange(n_sites)
    else:
        candidate = np.where(pos_mask)[0]

    # sort candidates by weight
    if direction == "high":
        order = np.argsort(site_weights[candidate])[::-1]
    elif direction == "low":
        order = np.argsort(site_weights[candidate])
    else:
        raise ValueError("direction must be 'high' or 'low'")

    selected = []
    cum_w = 0.0
    for j in order:
        idx = candidate[j]
        w = site_weights[idx]
        if w <= 0.0 and not include_zero_weights:
            continue
        selected.append(idx)
        cum_w += max(w, 0.0)
        if cum_w >= target_weight:
            break

    if not selected:
        raise ValueError("No sites selected (check coverage/include_zero_weights/weights).")

    # aggregate
    agg = None
    if normalize == "weights":
        den = 0.0
        for idx in selected:
            w = max(site_weights[idx], 0.0)
            mat = dist_matrices[idx]
            if agg is None:
                agg = (w * mat).copy()
            else:
                agg += w * mat
            den += w
        if den == 0.0:
            raise ValueError("Normalization denominator is zero (all selected weights <= 0).")
        agg /= den
    elif normalize == "count":
        for idx in selected:
            mat = dist_matrices[idx]
            if agg is None:
                agg = mat.copy()
            else:
                agg += mat
        agg /= len(selected)
    else:
        raise ValueError("normalize must be 'weights' or 'count'")

    return agg, selected
