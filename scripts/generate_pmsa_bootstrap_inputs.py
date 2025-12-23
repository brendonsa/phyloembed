#!/usr/bin/env python3
import argparse
from pathlib import Path
import random

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq


def scramble_alignment_columns(msa, frac=0.05, seed=None):
    """Scramble ~frac of columns in an aligned MSA (column-wise resampling with replacement)."""
    if seed is not None:
        random.seed(seed)

    n_seqs = len(msa)
    aln_len = msa.get_alignment_length()

    if frac == 1:
        positions = list(range(aln_len))
    else:
        n_cols = max(1, int(round(frac * aln_len)))
        n_cols = min(n_cols, aln_len)
        positions = random.sample(range(aln_len), n_cols)

    matrix = [list(str(rec.seq)) for rec in msa]

    for pos in positions:
        column = [matrix[i][pos] for i in range(n_seqs)]
        new_column = [random.choice(column) for _ in range(n_seqs)]
        for i in range(n_seqs):
            matrix[i][pos] = new_column[i]

    new_records = []
    for rec, row in zip(msa, matrix):
        seq_str = "".join(row)
        new_records.append(
            rec.__class__(id=rec.id, seq=Seq(seq_str), description="")
        )
    return new_records


def gap_postprocess(records, mode):
    """
    mode = "replace": replace gaps ('-') with 'X'
    mode = "remove":  delete gaps entirely
    """
    processed = []
    for rec in records:
        s = str(rec.seq)
        if mode == "replace":
            s2 = s.replace("-", "X")
        elif mode == "remove":
            s2 = s.replace("-", "")
        else:
            s2 = s

        processed.append(
            rec.__class__(id=rec.id, seq=Seq(s2), description="")
        )
    return processed


def build_pmsa_bootstrap(aligned_path, output_base, num, frac=0.05, suffix=""):
    msa = AlignIO.read(aligned_path, "fasta")

    base = Path(output_base)
    out_w = base / f"pmsabootstrapw{suffix}"
    out_wo = base / f"pmsabootstrapwithout{suffix}"
    out_w.mkdir(parents=True, exist_ok=True)
    out_wo.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Generating {num} PMSA bootstraps (w / without) → {output_base}{suffix}")
    for i in range(1, num + 1):
        scrambled = scramble_alignment_columns(msa, frac=frac, seed=1000 + i)

        # w: replace gaps with X
        recs_w = gap_postprocess(scrambled, "replace")
        rep_dir_w = out_w / str(i)
        rep_dir_w.mkdir(parents=True, exist_ok=True)
        out_path_w = rep_dir_w / "boot.fasta"
        SeqIO.write(recs_w, out_path_w, "fasta")
        print(f"[{i:03d}] pmsabootstrapw{suffix}  → {out_path_w}")

        # without: remove gaps altogether
        recs_wo = gap_postprocess(scrambled, "remove")
        rep_dir_wo = out_wo / str(i)
        rep_dir_wo.mkdir(parents=True, exist_ok=True)
        out_path_wo = rep_dir_wo / "boot.fasta"
        SeqIO.write(recs_wo, out_path_wo, "fasta")
        print(f"[{i:03d}] pmsabootstrapwithout{suffix} → {out_path_wo}")

    print("[DONE] PMSA bootstraps (w / without) generated.")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate PMSA-style bootstrap inputs: "
            "scramble a fraction of alignment columns and then either "
            "replace gaps with X (pmsabootstrapw) or remove gaps (pmsabootstrapwithout)."
        )
    )
    parser.add_argument(
        "-a", "--aligned", required=True, help="Input MSA FASTA file (aligned)"
    )
    parser.add_argument(
        "-o", "--output-base", required=True, help="Base output directory"
    )
    parser.add_argument(
        "-n", "--num", required=True, type=int, help="Number of bootstrap replicates"
    )
    parser.add_argument(
        "-f",
        "--fraction",
        type=float,
        default=0.05,
        help="Fraction of columns to scramble (default: 0.05)",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        default="",
        help="Suffix appended to output dirs (e.g. '_frac50')",
    )
    args = parser.parse_args()

    build_pmsa_bootstrap(
        args.aligned,
        args.output_base,
        args.num,
        frac=args.fraction,
        suffix=args.suffix,
    )


if __name__ == "__main__":
    main()
