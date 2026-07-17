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
    pooling: str = "mean",  # "mean", "min", "max", "concat", "cls"
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
            pooled = emb                      # (L, D)
        else:
            raise ValueError(f"Unknown pooling mode: {pooling}")

    # normalize for any non-concat pooling
    if pooling != "concat":
        norm = pooled.norm(p=2)
        if norm > 0:
            pooled = pooled / norm

    return pooled.cpu().numpy()


def process(input_fasta: str, output_path: str, gpu_id: int, pooling: str = "mean"):
    """
    If pooling != 'concat':
      FASTA -> pooled embeddings (N x D) -> CSV

    If pooling == 'concat':
      FASTA -> per-residue embeddings (N x L x D) -> NPZ (labels + embeddings)
      Requires all sequences have the same length (e.g., aligned FASTA).
    """
    pooling = pooling.lower()

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    client = ESMC.from_pretrained("esmc_600m").to(device)
    client.eval()

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    labels = []
    embs = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        labels.append(record.id)
        embs.append(embed_sequence(str(record.seq), client, device, pooling=pooling))

    # deterministic order
    order = np.argsort(np.asarray(labels, dtype=object))
    labels = np.asarray(labels, dtype=object)[order]
    embs = [embs[i] for i in order]

    if pooling != "concat":
        # 2D -> CSV
        df = pd.DataFrame(embs, index=labels)
        df.to_csv(
            output_file,
            index_label="strain",
            quoting=csv.QUOTE_ALL,
        )
        print(f"Saved {len(df)} embeddings → {output_file}")
        return

    # concat -> 3D -> NPZ
    try:
        embs_3d = np.stack(embs, axis=0)  # (N, L, D)
    except ValueError as e:
        raise ValueError(
            "pooling='concat' needs all sequences to have the same length to stack into (N, L, D). "
            "Use an aligned FASTA (same length for all sequences)."
        ) from e

    # enforce .npz extension
    if output_file.suffix.lower() != ".npz":
        output_file = output_file.with_suffix(".npz")

    np.savez_compressed(output_file, labels=labels, embeddings=embs_3d)
    print(f"Saved {len(labels)} raw embeddings → {output_file} (shape={embs_3d.shape})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute ESM-C embeddings.")
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV/NPZ path")
    parser.add_argument("--gpu-id", type=int, default=0, help="GPU device index (default: 0)")
    parser.add_argument(
        "--pooling",
        type=str,
        default="mean",
        choices=["mean", "min", "max", "concat", "cls", "avg", "average"],
        help="Pooling mode (default: mean). If concat, saves NPZ (N,L,D).",
    )
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id, pooling=args.pooling)