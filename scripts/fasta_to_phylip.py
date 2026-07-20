#!/usr/bin/env python3
"""Convert an aligned FASTA to sequential (single-row) relaxed PHYLIP."""
import argparse
from pathlib import Path


def read_fasta(path):
    name, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name = line[1:].strip().split()[0]
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description="Aligned FASTA -> sequential relaxed PHYLIP.")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    records = list(read_fasta(args.input))
    if not records:
        raise ValueError(f"No sequences read from {args.input}")

    lengths = {len(s) for _, s in records}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not equal length (found lengths {sorted(lengths)}); alignment required.")

    L = lengths.pop()
    with open(args.output, "w") as f:
        f.write(f"{len(records)} {L}\n")
        for name, seq in records:
            f.write(f"{name} {seq.upper()}\n")

    print(f"Wrote {len(records)} sequences (L={L}) -> {args.output}")


if __name__ == "__main__":
    main()
