import argparse
import pathlib
import re
from Bio import SeqIO
import torch
import numpy as np
import pandas as pd
from transformers import BertModel, BertTokenizer
import csv


def embed_residues(
    seq: str,
    tokenizer,
    model,
    device,
    keep_gaps: bool = False,
) -> torch.Tensor:
    # Clean sequence
    seq = seq.replace("*", "X")
    seq = re.sub(r"[UZOB]", "X", seq)
    if not keep_gaps:
        seq = seq.replace("-", "")

    if len(seq) == 0:
        raise ValueError("Sequence is empty after cleaning")

    tokens = tokenizer(
        " ".join(list(seq)),
        add_special_tokens=True,
        return_tensors="pt",
    ).to(device)

    with torch.no_grad():
        outputs = model(**tokens)
        emb = outputs.last_hidden_state[0, 1:len(seq) + 1]  # (L, D)

    return emb.cpu()


def embed_sequence(
    seq: str,
    tokenizer,
    model,
    device,
) -> torch.Tensor:
    return embed_residues(seq, tokenizer, model, device).mean(dim=0)


def process(input_fasta: str, output_path: str, gpu_id: int, pooling: str = "mean"):
    """Process FASTA → mean embeddings (L2-normalized) → CSV, or per-column → NPZ."""
    pooling = pooling.lower()

    # device management
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    # load model + tokenizer
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
    model = BertModel.from_pretrained("Rostlab/prot_bert").to(device)
    model.eval()

    output_file = pathlib.Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if pooling == "concat":
        labels, mats = [], []
        for record in SeqIO.parse(input_fasta, "fasta"):
            raw = str(record.seq)
            emb = embed_residues(raw, tokenizer, model, device, keep_gaps=True).numpy()
            if emb.shape[0] != len(raw):
                raise RuntimeError(
                    f"{record.id}: got {emb.shape[0]} tokens for {len(raw)} alignment columns."
                )
            labels.append(record.id)
            mats.append(emb.astype(np.float32))

        order = np.argsort(np.asarray(labels, dtype=object))
        labels = np.asarray(labels, dtype=object)[order]
        mats = [mats[i] for i in order]

        try:
            embs_3d = np.stack(mats, axis=0)
        except ValueError as e:
            raise ValueError(
                "pooling='concat' needs all sequences to have the same alignment length. Use an aligned FASTA."
            ) from e

        if output_file.suffix.lower() != ".npz":
            output_file = output_file.with_suffix(".npz")
        np.savez_compressed(output_file, labels=labels, embeddings=embs_3d)
        print(f"Saved {len(labels)} raw embeddings → {output_file} (shape={embs_3d.shape})")
        return

    # embed each record
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        label = record.id
        emb = embed_sequence(str(record.seq), tokenizer, model, device)  # torch.Tensor
        norm = emb.norm(p=2)
        if norm > 0:
            emb = emb / norm
        records.append((label, emb.numpy()))

    # convert to sorted dataframe
    records.sort(key=lambda x: x[0])
    df = pd.DataFrame([r[1] for r in records], index=[r[0] for r in records])

    df.to_csv(
        output_file,
        index_label="strain",
        quoting=csv.QUOTE_ALL,
    )

    print(f"Saved {len(df)} embeddings → {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute ProtBert embeddings (mean-pooled, L2-normalized)."
    )
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path")
    parser.add_argument("--gpu-id", type=int, default=0, help="GPU device index (default: 0)")
    parser.add_argument("--pooling", type=str, default="mean", choices=["mean", "concat"])
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id, pooling=args.pooling)
