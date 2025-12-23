import argparse
import pathlib
import re
from Bio import SeqIO
import torch
import pandas as pd
import csv
import numpy as np
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig


def embed_sequence(
    seq: str,
    client: ESMC,
    device,
    pooling: str = "cls",  # "mean", "min", "max", "concat", "cls"
) -> torch.Tensor:
    # seq = seq.replace("-", "")

    if len(seq) == 0:
        raise ValueError("Sequence is empty after cleaning")

    # run model
    protein = ESMProtein(sequence=seq)
    protein_tensor = client.encode(protein)

    with torch.no_grad():
        logits_output = client.logits(
            protein_tensor,
            LogitsConfig(sequence=True, return_embeddings=True),
        )

    emb = logits_output.embeddings
    if not isinstance(emb, torch.Tensor):
        emb = torch.tensor(emb)

    emb = emb.to(device).squeeze(0)  # (L, D)
    L, D = emb.shape
    seqlen = len(seq)

    # CLS pooling uses raw first token before stripping
    if pooling.lower() == "cls":
        if L in (seqlen + 1, seqlen + 2):
            pooled = emb[0]  # (D,)
        else:
            raise ValueError(
                f"Cannot safely identify CLS token: got L={L}, seq_len={seqlen}"
            )
    else:
        # strip special tokens for token-wise poolings
        if L == seqlen + 2:
            emb = emb[1:-1, :]  # [CLS] + tokens + [EOS]
        elif L == seqlen + 1:
            emb = emb[1:, :]    # [CLS] + tokens

        if pooling in ("mean", "avg", "average"):
            pooled = emb.mean(dim=0)          # (D,)
        elif pooling == "max":
            pooled, _ = emb.max(dim=0)        # (D,)
        elif pooling == "min":
            pooled, _ = emb.min(dim=0)        # (D,)
        elif pooling == "concat":
            pooled = emb                      # (L, D)  <-- per-residue embeddings
        else:
            raise ValueError(f"Unknown pooling mode: {pooling}")

    # normalize for any non-concat pooling
    if pooling != "concat":
        norm = pooled.norm(p=2)
        if norm > 0:
            pooled = pooled / norm

    return pooled.cpu().numpy()


def process(input_fasta: str, output_path: str, gpu_id: int):
    """Process FASTA → mean embeddings → CSV (ESM-C backend, same I/O as original)."""
    # device management
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    # load ESM-C model client
    client = ESMC.from_pretrained("esmc_600m").to(device)
    client.eval()

    # embed each record
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        label = record.id
        emb = embed_sequence(str(record.seq), client, device)
        records.append((label, emb))

    # convert to sorted dataframe (same layout as original)
    records.sort(key=lambda x: x[0])
    df = pd.DataFrame([r[1] for r in records], index=[r[0] for r in records])

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        output_file,
        index_label="strain",
        quoting=csv.QUOTE_ALL,  # ensure names with spaces are preserved
    )

    print(f"Saved {len(df)} embeddings → {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute ESM-C embeddings (chunk-aware, 2048 context)."
    )
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path")
    parser.add_argument("--gpu-id", type=int, default=0, help="GPU device index (default: 0)")
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id)
