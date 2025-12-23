import argparse
import numpy as np
import pandas as pd

def pairwise_distances(X: np.ndarray, metric: str) -> np.ndarray:
    if metric == "cosine":
        Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
        D = 1.0 - np.dot(Xn, Xn.T)
    elif metric == "euclidean":
        r = np.sum(X * X, axis=1, keepdims=True)
        D = np.sqrt(np.maximum(r + r.T - 2.0 * np.dot(X, X.T), 0.0))
    else:
        raise ValueError(f"Unsupported metric: {metric}")
    np.fill_diagonal(D, 0.0)
    return D

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compute pairwise distance matrix.")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-m", "--metric", required=True, choices=["cosine", "euclidean"])
    args = p.parse_args()

    df = pd.read_csv(args.input, index_col=0)
    X = df.to_numpy(dtype=np.float32, copy=False)

    D = pairwise_distances(X, args.metric)
    pd.DataFrame(D, index=df.index, columns=df.index).to_csv(args.output)
    print(f"Saved {args.metric} distances {args.output}")
