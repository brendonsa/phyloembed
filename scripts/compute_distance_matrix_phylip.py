import argparse
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances

def pairwise_distances(
    X: np.ndarray,
    metric: str,
    tiny_tol: float = 1e-10,
    eps: float = 1e-8,
) -> np.ndarray:
    if metric == "cosine":
        D = cosine_distances(X)
        # Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
        # D = 1.0 - np.dot(Xn, Xn.T)

    elif metric == "euclidean":
       D = euclidean_distances(X)

        # r = np.sum(X * X, axis=1, keepdims=True)
        # D = np.sqrt(np.maximum(r + r.T - 2.0 * np.dot(X, X.T), 0.0))

    elif metric == "correlation":
        # center each vector by its own mean
        Xc = X - X.mean(axis=1, keepdims=True)
        Xc_norm = Xc / (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)
        # correlation distance = 1 - corr
        D = 1.0 - np.dot(Xc_norm, Xc_norm.T)

    elif metric == "mahalanobis":
        # global centering
        Xc = X - X.mean(axis=0, keepdims=True)

        # covariance over features
        cov = np.cov(Xc, rowvar=False)

        # small Tikhonov regularization for stability
        d = cov.shape[0]
        reg = 1e-6 * (np.trace(cov) / d)
        cov = cov + reg * np.eye(d, dtype=cov.dtype)

        # whitening transform using eigendecomposition
        eigvals, eigvecs = np.linalg.eigh(cov)
        inv_sqrt = 1.0 / np.sqrt(np.maximum(eigvals, 1e-12))
        # W = Q * diag(lambda^{-1/2})
        W = eigvecs * inv_sqrt

        # whitened features
        Xw = Xc @ W

        # Mahalanobis distance = Euclidean distance in whitened space
        r = np.sum(Xw * Xw, axis=1, keepdims=True)
        D = np.sqrt(np.maximum(r + r.T - 2.0 * np.dot(Xw, Xw.T), 0.0))

    else:
        raise ValueError(f"Unsupported metric: {metric}")

    # clamp negatives / -0.0 from floating point
    np.maximum(D, 0.0, out=D)

    # enforce exact zeros on the diagonal
    np.fill_diagonal(D, 0.0)

    # replace negligible *off-diagonal* distances with epsilon
    n = D.shape[0]
    off_diag = ~np.eye(n, dtype=bool)
    tiny = (D < tiny_tol) & off_diag
    D[tiny] = eps

    return D



def write_phylip_relaxed(D: np.ndarray, labels, path: str, precision: int = 6) -> None:
    n = len(labels)
    with open(path, "w") as f:
        # first line: number of taxa
        f.write(f"{n}\n")
        # relaxed PHYLIP: full names, whitespace-separated
        for i, name in enumerate(labels):
            parts = [str(name)]
            parts.extend(f"{D[i, j]:.{precision}f}" for j in range(n))
            f.write(" ".join(parts) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compute pairwise distance matrix in relaxed PHYLIP format.")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-m", "--metric", required=True, choices=["cosine", "euclidean", "mahalanobis", "correlation"])
    args = p.parse_args()

    df = pd.read_csv(args.input, index_col=0)
    X = df.to_numpy(dtype=np.float32, copy=False)

    D = pairwise_distances(X, args.metric)
    labels = df.index.tolist()

    write_phylip_relaxed(D, labels, args.output)
    print(f"Saved {args.metric} distances in relaxed PHYLIP format to {args.output}")
