#!/usr/bin/env python3
import argparse
from pathlib import Path
import random
from Bio import AlignIO, SeqIO, Phylo
from Bio.Phylo.Consensus import bootstrap, get_support
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def mutate_unaligned(records, mutation_rate=0.01, seed=None):
    """Generate pseudo-bootstrap replicate by random residue mutations in unaligned sequences."""
    if seed is not None:
        random.seed(seed)

    mutated_records = []
    for rec in records:
        seq_list = list(str(rec.seq))
        n_mut = max(1, int(len(seq_list) * mutation_rate))
        positions = random.sample(range(len(seq_list)), n_mut)
        for pos in positions:
            aa = seq_list[pos]
            if aa.isalpha():
                new_aa = aa
                while new_aa == aa:
                    new_aa = random.choice("ACDEFGHIKLMNPQRSTVWY")
                seq_list[pos] = new_aa
        new_rec = rec.__class__(id=rec.id, seq="".join(seq_list), description="")
        mutated_records.append(new_rec)
    return mutated_records


def build_bootstrap(msa_path, output_dir, num, metric):
    """Classic MSA bootstrap."""
    msa = AlignIO.read(msa_path, "fasta")
    base = Path(output_dir)
    base.mkdir(parents=True, exist_ok=True)
    calculator = DistanceCalculator(metric)
    constructor = DistanceTreeConstructor(calculator, "nj")

    print(f"[INFO] Generating {num} classical bootstraps → {output_dir}")
    trees = []
    for i, boot_alignment in enumerate(bootstrap(msa, num), start=1):
        rep_dir = base / str(i)
        rep_dir.mkdir(parents=True, exist_ok=True)
        out_path = rep_dir / "boot.fasta"    
        SeqIO.write(boot_alignment, out_path, "fasta")
        print(f"[{i:03d}] Saved bootstrap MSA → {out_path}")

    # base_tree = constructor.build_tree(msa)
    # tree_path = base / "tree_nj.nwk"
    # Phylo.write(base_tree, tree_path, "newick")
    # print(f"[DONE] Bootstrap NJ tree saved → {tree_path}")


def build_pseudo_bootstrap(fasta_path, output_dir, num, mutation_rate=0.01):
    """Pseudo-bootstraps for unaligned sequences."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    base = Path(output_dir)
    base.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Generating {num} pseudo-bootstraps (unaligned) → {output_dir}")
    random.seed(42)
    for i in range(1, num + 1):
        rep_dir = base / str(i)
        rep_dir.mkdir(parents=True, exist_ok=True)
        out_path = rep_dir / "boot.fasta"
        mutated_records = mutate_unaligned(records, mutation_rate=mutation_rate, seed=42 + i)
        SeqIO.write(mutated_records, out_path, "fasta")
        print(f"[{i:03d}] Saved pseudo-bootstrap FASTA → {out_path}")

    print(f"[DONE] Generated {num} pseudo-bootstrap replicates (unaligned).")


def main(msa_input, unaligned_input, output_base, num, metric):
    msa_out = Path(output_base) / "bootstrap"
    pseudo_out = Path(output_base) / "pseudo_bootstrap"

    build_bootstrap(msa_input, msa_out, num, metric)
    build_pseudo_bootstrap(unaligned_input, pseudo_out, num)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate both classical (MSA-based) and pseudo-bootstrap (unaligned) replicates."
    )
    parser.add_argument("-a", "--aligned", required=True, help="Input MSA FASTA file (aligned)")
    parser.add_argument("-u", "--unaligned", required=True, help="Input FASTA file (unaligned)")
    parser.add_argument("-o", "--output-base", required=True, help="Base output directory")
    parser.add_argument("-n", "--num", required=True, type=int, help="Number of bootstrap replicates")
    parser.add_argument(
        "-m",
        "--metric",
        default="blosum62",
        choices=[
            "identity",
            "blastn",
            "trans",
            "benner6",
            "benner22",
            "benner74",
            "grantham",
            "pam250",
            "blosum62",
        ],
        help="Distance metric for NJ construction (default: blosum62)",
    )
    args = parser.parse_args()
    main(args.aligned, args.unaligned, args.output_base, args.num, args.metric)
