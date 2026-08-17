import argparse
import csv
import pathlib

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
import esm


def load_msa(input_fasta):
    labels, seqs = [], []
    for record in SeqIO.parse(input_fasta, "fasta"):
        labels.append(record.id)
        seqs.append(str(record.seq))
    return labels, seqs


def process(input_fasta, output_path, gpu_id, pooling="mean"):
    pooling = pooling.lower()

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    labels, seqs = load_msa(input_fasta)
    data = [(labels[i], seqs[i]) for i in range(len(labels))]

    _, _, batch_tokens = batch_converter([data])  # (1, N, L+1)
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        out = model(batch_tokens, repr_layers=[12])
    reps = out["representations"][12][0]  # (N, L+1, D)
    reps = reps[:, 1:, :]  # drop leading BOS token -> (N, L, D)

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    order = np.argsort(np.asarray(labels, dtype=object))
    labels_sorted = np.asarray(labels, dtype=object)[order]
    reps = reps[order]

    if pooling != "concat":
        pooled = reps.mean(dim=1)  # (N, D)
        norm = pooled.norm(p=2, dim=1, keepdim=True)
        pooled = pooled / norm.clamp_min(1e-12)
        df = pd.DataFrame(pooled.cpu().numpy(), index=labels_sorted)
        df.to_csv(output_file, index_label="strain", quoting=csv.QUOTE_ALL)
        print(f"Saved {len(df)} embeddings -> {output_file}")
        return

    if output_file.suffix.lower() != ".npz":
        output_file = output_file.with_suffix(".npz")
    np.savez_compressed(output_file, labels=labels_sorted, embeddings=reps.cpu().numpy())
    print(f"Saved {len(labels_sorted)} raw embeddings -> {output_file} (shape={tuple(reps.shape)})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute MSA Transformer embeddings.")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--pooling", type=str, default="mean", choices=["mean", "concat"])
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id, pooling=args.pooling)
