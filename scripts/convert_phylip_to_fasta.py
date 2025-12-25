#!/usr/bin/env python3
import argparse
from pathlib import Path

from Bio import AlignIO, SeqIO


def read_alignment(path: str):
    try:
        return AlignIO.read(path, "phylip-relaxed")
    except Exception:
        return AlignIO.read(path, "phylip")


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert PHYLIP alignment to FASTA.")
    ap.add_argument("-i", "--input", required=True, help="Input PHYLIP alignment.")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA alignment.")
    args = ap.parse_args()

    aln = read_alignment(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(aln, out_path, "fasta")


if __name__ == "__main__":
    main()
