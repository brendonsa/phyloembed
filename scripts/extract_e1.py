import argparse
import pathlib
import csv
import numpy as np

from Bio import SeqIO
import torch
import pandas as pd

from E1.modeling import E1ForMaskedLM
from E1.predictor import E1Predictor


def clean_sequence(seq: str) -> str:
    seq = seq.replace("-", "")
    return seq


def embed_e1(
    seq: str,
    predictor: E1Predictor,
    seq_id: str,
    pooling: str = "mean",  # "mean" or "concat"
) -> np.ndarray:
    """
    pooling='mean'   -> returns (D,)
    pooling='concat' -> returns (L,D)
    """
    pooling = pooling.lower()
    if pooling not in ("mean", "concat"):
        raise ValueError("pooling must be 'mean' or 'concat'")

    # run predictor on a single sequence (keep same API style as ESM-C file)
    pred_iter = predictor.predict(
        sequences=[seq],
        sequence_ids=[seq_id],
        context_seqs=None,
    )
    pred = next(iter(pred_iter))

    if pooling == "mean":
        emb = pred["mean_token_embeddings"]  # (D,)
        norm = emb.norm(p=2)
        if norm > 0:
            emb = emb / norm
        return emb.detach().cpu().numpy()

    # concat
    emb = pred["token_embeddings"]  # (L, D)
    return emb.detach().cpu().numpy()


def process(input_fasta: str, output_path: str, gpu_id: int, pooling: str = "mean", max_batch_tokens: int = 16384):
    """
    If pooling != 'concat':
      FASTA -> pooled embeddings (N x D) -> CSV

    If pooling == 'concat':
      FASTA -> per-token embeddings (N x L x D) -> NPZ (labels + embeddings)
      Requires all sequences same length to stack.
    """
    pooling = pooling.lower()

    # device
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    model_name = "Profluent-Bio/E1-600m"

    # model
    model = E1ForMaskedLM.from_pretrained(model_name, torch_dtype=torch.float32).to(device).eval()

    # predictor: request the right fields based on pooling
    fields_to_save = ["token_embeddings"] if pooling == "concat" else ["mean_token_embeddings"]
    predictor = E1Predictor(
        model=model,
        max_batch_tokens=max_batch_tokens,
        fields_to_save=fields_to_save,
    )

    # embed each record
    labels = []
    embs = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = clean_sequence(str(record.seq))
        if len(seq) == 0:
            print(f"Skipping {record.id}: empty after cleaning")
            continue

        labels.append(record.id)
        embs.append(embed_e1(seq, predictor, record.id, pooling=pooling))

    print(f"Embedded {len(embs)} sequences from {input_fasta}")

    # deterministic order
    order = np.argsort(np.asarray(labels, dtype=object))
    labels = np.asarray(labels, dtype=object)[order]
    embs = [embs[i] for i in order]

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if pooling != "concat":
        # (N, D) -> CSV
        df = pd.DataFrame(embs, index=labels)
        df.to_csv(
            output_file,
            index_label="strain",
            quoting=csv.QUOTE_ALL,
        )
        print(f"Saved {len(df)} embeddings → {output_file}")
        return

    # (N, L, D) -> NPZ
    try:
        embs_3d = np.stack(embs, axis=0)
    except ValueError as e:
        raise ValueError(
            "pooling='concat' needs all sequences to have the same length to stack into (N, L, D). "
            "Use an aligned FASTA (same length for all sequences)."
        ) from e

    if output_file.suffix.lower() != ".npz":
        output_file = output_file.with_suffix(".npz")

    np.savez_compressed(output_file, labels=labels, embeddings=embs_3d)
    print(f"Saved {len(labels)} raw embeddings → {output_file} (shape={embs_3d.shape})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute E1 embeddings (single-sequence mode)."
    )
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV/NPZ path")
    parser.add_argument("--gpu-id", type=int, default=0, help="GPU device index (default: 0)")
    parser.add_argument(
        "--pooling",
        type=str,
        default="mean",
        choices=["mean", "concat"],
        help="Pooling mode (default: mean). If concat, saves NPZ (N,L,D).",
    )
    parser.add_argument(
        "--max-batch-tokens",
        type=int,
        default=4096,
        help="Batch token budget for E1Predictor (default: 512)",
    )
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id, pooling=args.pooling, max_batch_tokens=args.max_batch_tokens)
