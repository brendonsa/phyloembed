import argparse
import csv
import pathlib

import numpy as np
import pandas as pd
from Bio import SeqIO


def load_ngram_vectors(path):
    table = pd.read_csv(path, sep="\t", index_col=0)
    vectors = {ngram: row.to_numpy(dtype=np.float64) for ngram, row in table.iterrows()}
    return vectors, vectors["<unk>"]


def embed_sequence(seq, vectors, unk, n=3):
    """Whole-sequence ProtVec: sum of 3-gram vectors over all n reading frames."""
    seq = seq.upper().replace("-", "").replace(".", "")
    total = np.zeros_like(unk)
    for offset in range(n):
        for i in range(offset, len(seq) - n + 1, n):
            total += vectors.get(seq[i:i + n], unk)
    return total


def embed_columns(seq, vectors, unk, n=3):
    """One vector per alignment column: the overlapping n-gram starting at that column."""
    seq = seq.upper()
    L = len(seq)
    return np.stack([
        vectors.get(seq[i:i + n], unk) if i + n <= L else unk
        for i in range(L)
    ])


def process(input_fasta, vectors_path, output_path, pooling="mean"):
    pooling = pooling.lower()
    vectors, unk = load_ngram_vectors(vectors_path)

    labels, embs = [], []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if pooling == "concat":
            embs.append(embed_columns(str(record.seq), vectors, unk))
        else:
            embs.append(embed_sequence(str(record.seq), vectors, unk))
        labels.append(record.id)

    order = np.argsort(np.asarray(labels, dtype=object))
    labels = np.asarray(labels, dtype=object)[order]
    embs = [embs[i] for i in order]

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if pooling != "concat":
        pooled = np.stack(embs, axis=0)
        norms = np.linalg.norm(pooled, axis=1, keepdims=True)
        pooled = pooled / np.clip(norms, 1e-12, None)
        df = pd.DataFrame(pooled, index=labels)
        df.to_csv(output_file, index_label="strain", quoting=csv.QUOTE_ALL)
        print(f"Saved {len(df)} embeddings -> {output_file}")
        return

    try:
        embs_3d = np.stack(embs, axis=0)
    except ValueError as e:
        raise ValueError(
            "pooling='concat' needs all sequences to have the same alignment length. Use an aligned FASTA."
        ) from e

    if output_file.suffix.lower() != ".npz":
        output_file = output_file.with_suffix(".npz")
    np.savez_compressed(output_file, labels=labels, embeddings=embs_3d)
    print(f"Saved {len(labels)} raw embeddings -> {output_file} (shape={embs_3d.shape})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute ProtVec embeddings from pretrained 3-gram vectors.")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--vectors", required=True, help="protVec_100d_3grams.csv (tab-separated)")
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--pooling", type=str, default="mean", choices=["mean", "concat"])
    args = parser.parse_args()

    process(args.input, args.vectors, args.output, pooling=args.pooling)
