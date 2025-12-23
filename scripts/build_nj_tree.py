#!/usr/bin/env python3
import argparse
import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import nj
from skbio.io import write

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def build_nj(input_csv: str, output_newick: str):
    df = pd.read_csv(input_csv, index_col=0)
    names = list(df.index)
    matrix = df.values.astype(float)
    dm = DistanceMatrix(matrix, ids=names)
    tree = nj(dm)
    with open(output_newick, "w") as f:
        write(tree, format="newick", into=f)
    print(f"NJ tree written to {output_newick}")

def build_nj_from_msa(input_msa: str, output_newick: str, metric: str):
    msa = AlignIO.read(input_msa, "fasta")
    calculator = DistanceCalculator(metric)
    constructor = DistanceTreeConstructor(calculator, "nj")
    dm = calculator.get_distance(msa)
    tree = constructor.nj(dm)
    Phylo.write(tree, output_newick, "newick")
    print(f"NJ tree built from MSA → {output_newick}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--from-msa",
        action="store_true",
        help="Treat input as FASTA MSA instead of distance matrix CSV",
    )
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
        help="Distance metric for MSA-based NJ (default: blosum62)",
    )
    args = parser.parse_args()

    if args.from_msa:
        build_nj_from_msa(args.input, args.output, args.metric)
    else:
        build_nj(args.input, args.output)
