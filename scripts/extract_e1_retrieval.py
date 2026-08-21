import argparse
import csv
import pathlib
from collections import Counter

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO

from E1.modeling import E1ForMaskedLM
from E1.predictor import E1Predictor


def kmer_profile(seqs, k=3):
    vocab = {}
    counters = []
    for s in seqs:
        c = Counter(s[j:j + k] for j in range(len(s) - k + 1))
        counters.append(c)
        for kmer in c:
            vocab.setdefault(kmer, len(vocab))
    mat = np.zeros((len(seqs), len(vocab)), dtype=np.float32)
    for i, c in enumerate(counters):
        for kmer, cnt in c.items():
            mat[i, vocab[kmer]] = cnt
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return mat / norms


def build_context(clean_seqs, i, budget, sim_matrix):
    """Build context only for the budget (in token size).
       This is used to reduce the memory load for OOM issue.
       First redundant(most similar) sequences are dropped"""
    idxs = np.array([j for j in range(len(clean_seqs)) if j != i and len(clean_seqs[j]) > 0])
    lengths = np.array([len(clean_seqs[j]) for j in idxs])
    total = int(lengths.sum())
    if total <= budget:
        return ",".join(clean_seqs[j] for j in idxs)

    sub = sim_matrix[np.ix_(idxs, idxs)].copy()
    np.fill_diagonal(sub, -np.inf)
    kept = np.ones(len(idxs), dtype=bool)

    while total > budget:
        a, b = np.unravel_index(np.argmax(sub), sub.shape)
        drop = a if lengths[a] >= lengths[b] else b
        kept[drop] = False
        total -= lengths[drop]
        sub[drop, :] = -np.inf
        sub[:, drop] = -np.inf

    return ",".join(clean_seqs[j] for j in idxs[kept])


def embed_e1_retrieval(seq, context, predictor, seq_id, pooling="mean"):
    """
    pooling='mean'   -> returns (D,)
    pooling='concat' -> returns (L,D), L = len(seq)
    """
    context_seqs = {"family": context} if context else None
    pred = next(iter(predictor.predict(
        sequences=[seq],
        sequence_ids=[seq_id],
        context_seqs=context_seqs,
    )))

    if pooling == "mean":
        emb = pred["mean_token_embeddings"]
        norm = emb.norm(p=2)
        if norm > 0:
            emb = emb / norm
        return emb.detach().cpu().numpy()

    emb = pred["token_embeddings"]  # (L_total, D) -- context+query if context given
    emb = emb[-len(seq):]  # keep only the query's own positions
    if emb.shape[0] != len(seq):
        raise RuntimeError(
            f"{seq_id}: expected {len(seq)} query positions, got {emb.shape[0]} "
            "-- context/query slicing assumption may be wrong for this model version."
        )
    return emb.detach().cpu().numpy()


def process(input_fasta, output_path, gpu_id, pooling="mean", max_batch_tokens=4096,
            oom_max_context=13000, oom_min_context=1000, oom_backoff=0.5):
    pooling = pooling.lower()

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    model_name = "Profluent-Bio/E1-600m"
    model = E1ForMaskedLM.from_pretrained(model_name, torch_dtype=torch.float32).to(device).eval()
    fields_to_save = ["token_embeddings"] if pooling == "concat" else ["mean_token_embeddings"]
    predictor = E1Predictor(
        model=model,
        max_batch_tokens=max_batch_tokens,
        fields_to_save=fields_to_save,
    )

    records = list(SeqIO.parse(input_fasta, "fasta"))
    ids = [r.id for r in records]
    raw_seqs = [str(r.seq) for r in records]
    clean_seqs = [s.replace("-", "") for s in raw_seqs]

    kmer_mat = kmer_profile(clean_seqs)
    sim_matrix = kmer_mat @ kmer_mat.T

    labels, embs = [], []
    for i, (seq_id, seq) in enumerate(zip(ids, clean_seqs)):
        if len(seq) == 0:
            print(f"Skipping {seq_id}: empty after removing gaps")
            continue

        budget = float("inf")  # first attempt: full family context, no dropping
        while True:
            context = build_context(clean_seqs, i, budget, sim_matrix)
            try:
                emb = embed_e1_retrieval(seq, context, predictor, seq_id, pooling=pooling)
                break
            except RuntimeError as e:
                if "out of memory" not in str(e).lower():
                    raise
                torch.cuda.empty_cache()
                budget = oom_max_context if budget == float("inf") else int(budget * oom_backoff)
                if budget < oom_min_context:
                    raise
                print(f"{seq_id}: OOM, retrying with context budget={budget}")

        if pooling == "concat":
            # re-expand to full alignment length; gaps get a tiny non-zero filler
            # (a true zero vector makes cosine distance divide by zero -> nan)
            full = np.full((len(raw_seqs[i]), emb.shape[1]), 1e-8, dtype=np.float32)
            pos = 0
            for k, ch in enumerate(raw_seqs[i]):
                if ch != "-":
                    full[k] = emb[pos]
                    pos += 1
            emb = full

        labels.append(seq_id)
        embs.append(emb)

    print(f"Embedded {len(embs)} sequences (retrieval-augmented, whole-family context) from {input_fasta}")

    order = np.argsort(np.asarray(labels, dtype=object))
    labels = np.asarray(labels, dtype=object)[order]
    embs = [embs[i] for i in order]

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if pooling != "concat":
        df = pd.DataFrame(embs, index=labels)
        df.to_csv(output_file, index_label="strain", quoting=csv.QUOTE_ALL)
        print(f"Saved {len(df)} embeddings -> {output_file}")
        return

    try:
        embs_3d = np.stack(embs, axis=0)
    except ValueError as e:
        raise ValueError(
            "pooling='concat' needs all sequences to have the same alignment length to stack into (N, L, D). "
            "Use an aligned FASTA."
        ) from e

    if output_file.suffix.lower() != ".npz":
        output_file = output_file.with_suffix(".npz")
    np.savez_compressed(output_file, labels=labels, embeddings=embs_3d)
    print(f"Saved {len(labels)} raw embeddings -> {output_file} (shape={embs_3d.shape})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute E1 embeddings in retrieval-augmented mode (whole-family context)."
    )
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--pooling", type=str, default="mean", choices=["mean", "concat"])
    parser.add_argument("--oom-max-context", type=int, default=13000)
    parser.add_argument("--oom-min-context", type=int, default=1000)
    parser.add_argument("--oom-backoff", type=float, default=0.5)
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id, pooling=args.pooling,
            oom_max_context=args.oom_max_context,
            oom_min_context=args.oom_min_context,
            oom_backoff=args.oom_backoff)
