import argparse
import pathlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

NUC_CHARS = set("ACGTUNacgtun")  # Typical nucleotide characters

def is_nucleotide_sequence(seq: str) -> bool:
    """Return True if sequence looks like nucleotide (heuristic)."""
    seq_set = set(seq)
    return seq_set.issubset(NUC_CHARS)

def translate_fasta(input_fasta: str, output_fasta: str):
    """Translate nucleotide sequences to amino acids if needed."""
    records_in = list(SeqIO.parse(input_fasta, "fasta"))

    if all(is_nucleotide_sequence(str(rec.seq)) for rec in records_in):
        print("Detected nucleotide sequences. Translating to amino acids...")
        
        records_out = []
        for record in records_in:
            aa_seq = record.seq.translate(to_stop=False)
            records_out.append(SeqRecord(aa_seq, id=record.id, description=""))

    else:
        print("Detected amino acid sequences. Copying without translation...")
        records_out = [SeqRecord(rec.seq, id=rec.id, description="") for rec in records_in]

    out_path = pathlib.Path(output_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records_out, out_path, "fasta")

    print(f"Saved {len(records_out)} sequences → {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate nucleotide FASTA to amino acid FASTA if needed.")
    parser.add_argument("--input", "-i", required=True, help="Input FASTA file (nucleotide or AA)")
    parser.add_argument("--output", "-o", required=True, help="Output amino acid FASTA file")
    args = parser.parse_args()

    translate_fasta(args.input, args.output)
