import argparse
import pathlib
import re
from Bio import SeqIO
import torch
import pandas as pd
from transformers import T5Tokenizer, T5EncoderModel


def embed_sequence(
    seq: str,
    tokenizer,
    model,
    device,
    max_len: int = 512,
) -> torch.Tensor:
    # Clean sequence
    seq = seq.replace("*", "X")
    seq = re.sub(r"[UZOB]", "X", seq)
    seq = seq.replace("-", "")  # gaps → X

    if len(seq) == 0:
        raise ValueError("Sequence is empty after cleaning")

    tokens = tokenizer(
        " ".join(list(seq)),
        add_special_tokens=True,
        return_tensors="pt",
        truncation=True,
        max_length=max_len,
    ).to(device)

    with torch.no_grad():
        outputs = model(**tokens)
        # per-residue embeddings (exclude first special token)
        emb = outputs.last_hidden_state[0, 1:len(seq) + 1]  # (L, D)
        pooled = emb.mean(dim=0)  # (D,)

    return pooled.cpu()  # torch tensor on CPU


def process(input_fasta: str, output_path: str, gpu_id: int):
    """Compute ProtT5 embeddings and save as DataFrame."""
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{gpu_id}")
        print(f"Using GPU: {torch.cuda.get_device_name(gpu_id)}")
    else:
        device = torch.device("cpu")
        print("CUDA not available, using CPU")

    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(device)
    model.eval()

    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        label = record.id
        emb = embed_sequence(str(record.seq), tokenizer, model, device)  # torch.Tensor
        norm = emb.norm(p=2)
        if norm > 0:
            emb = emb / norm
        records.append((label, emb.numpy()))

    records.sort(key=lambda x: x[0])
    df = pd.DataFrame([r[1] for r in records], index=[r[0] for r in records])

    out_path = pathlib.Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path)

    print(f"Saved {len(df)} ProtT5 embeddings → {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute ProtT5 embeddings (mean-pooled, L2-normalized)."
    )
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path")
    parser.add_argument("--gpu-id", type=int, default=0, help="GPU device index (default: 0)")
    args = parser.parse_args()

    process(args.input, args.output, args.gpu_id)
