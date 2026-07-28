#!/usr/bin/env python3
import argparse
import importlib.util
import os
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description="Infer a phylogenetic tree with Phyla.")
    ap.add_argument("--input", required=True, help="Protein FASTA (aligned or unaligned).")
    ap.add_argument("--output", required=True, help="Output Newick (.nwk) file.")
    ap.add_argument("--model", default="phyla-beta", help="Phyla model name (default: phyla-beta).")
    ap.add_argument("--device", default="cuda:0", help="Torch device (default: cuda:0).")
    args = ap.parse_args()

    # model.py does `from utils.utils import ...`, needs phyla/'s own dir on sys.path
    phyla_pkg_dir = os.path.dirname(importlib.util.find_spec("phyla").origin)
    sys.path.insert(0, phyla_pkg_dir)

    from phyla.utils.eval_configs import Config
    from phyla import phyla as PhylaModel

    config = Config()
    config.model.model_name = args.model
    config.model.n_layer = 16

    model = PhylaModel(config, device=args.device)
    ckpt_paths = {"phyla-beta": Path("weights/11564369"),
                  "phyla-alpha": Path("weights/phyla_alpha_291M_state_dict.pt")}
    ckpt_path = ckpt_paths.get(model.version)
    if ckpt_path is not None and ckpt_path.exists():
        model.load(checkpoint_file=str(ckpt_path))
    else:
        model.load()
    if args.device.startswith("cuda"):
        model = model.cuda()

    encoded_aa, cls_token_mask, sequence_mask, sequence_names = model.encode_fasta(args.input)
    preds = model(encoded_aa, sequence_mask, cls_token_mask)
    tree = model.reconstruct_tree(preds, sequence_names)

    tree.write(args.output)
    print(f"[phyla] wrote tree -> {args.output}")


if __name__ == "__main__":
    main()
