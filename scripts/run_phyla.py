#!/usr/bin/env python3
import argparse


def main():
    ap = argparse.ArgumentParser(description="Infer a phylogenetic tree with Phyla.")
    ap.add_argument("--input", required=True, help="Protein FASTA (aligned or unaligned).")
    ap.add_argument("--output", required=True, help="Output Newick (.nwk) file.")
    ap.add_argument("--model", default="phyla-beta", help="Phyla model name (default: phyla-beta).")
    ap.add_argument("--device", default="cuda:0", help="Torch device (default: cuda:0).")
    args = ap.parse_args()

    from phyla import phyla

    model = phyla(name=args.model, device=args.device).load()
    if args.device.startswith("cuda"):
        model = model.cuda()

    encoded_aa, cls_token_mask, sequence_mask, sequence_names = model.encode_fasta(args.input)
    preds = model(encoded_aa, sequence_mask, cls_token_mask)
    tree = model.reconstruct_tree(preds, sequence_names)

    tree.write(args.output)
    print(f"[phyla] wrote tree -> {args.output}")


if __name__ == "__main__":
    main()
