#!/usr/bin/env python3
"""Run NeuralNJ (https://github.com/ZhangXinru99/NeuralNJ) on a PHYLIP alignment and
write a Newick tree.

NeuralNJ inference (`finetune_rl_search.py --infer_opt Argmax`) scans a directory of
`.phy` files and writes one `<name>.tre` per input under
`output/{opt}_dim{embed_dim}_patch{patch_size}/{instance_dirname}/`.

NOTE: NeuralNJ only models nucleotides (`DNA_WITH_GAP`, vocab size 4). Feed it a
nucleotide alignment (e.g. `aligned_nucleotides.fasta`), not a protein MSA.

Assumes the NeuralNJ repo is checked out as a sibling directory (default: ./NeuralNJ),
its conda env is active, and raxmlpy / IQTree / RAxML are installed per its README.
"""
import argparse
import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

import yaml


def main():
    ap = argparse.ArgumentParser(description="Infer a phylogenetic tree with NeuralNJ.")
    ap.add_argument("--phy", required=True, help="Input alignment in PHYLIP (.phy) format.")
    ap.add_argument("--output", required=True, help="Output Newick (.nwk) file.")
    ap.add_argument("--repo", default="NeuralNJ", help="Path to the NeuralNJ checkout (default: NeuralNJ).")
    ap.add_argument("--config-template", default=None,
                    help="Config YAML to base inference on (default: <repo>/config/finetune_reinforce_search_example.yaml).")
    ap.add_argument("--checkpoint", default=None,
                    help="Override reload_checkpoint_path (relative to repo). Default: keep template value.")
    ap.add_argument("--infer-opt", default="Argmax", choices=["Argmax", "Search", "Finetune"])
    ap.add_argument("--evolution-model", default="GTR+I+G")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    template = Path(args.config_template) if args.config_template else repo / "config" / "finetune_reinforce_search_example.yaml"
    cfg = yaml.safe_load(Path(template).read_text())

    # Isolated instance directory (relative to repo, since NeuralNJ prefixes with its own dir).
    inst_name = f"_pe_instance_{uuid.uuid4().hex[:8]}"
    inst_dir = repo / inst_name
    inst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(args.phy, inst_dir / "aln.phy")

    cfg["instance_path"] = f"./{inst_name}"
    if args.checkpoint:
        cfg["reload_checkpoint_path"] = args.checkpoint

    embed_dim = cfg["model"]["embed_dim"]
    patch_size = cfg["model"]["patch_size"]

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tf:
        yaml.safe_dump(cfg, tf)
        cfg_path = tf.name

    try:
        subprocess.run(
            ["python", "finetune_rl_search.py",
             "--config", cfg_path,
             "--infer_opt", args.infer_opt,
             "--evolution_model", args.evolution_model],
            cwd=str(repo),
            check=True,
        )

        tre = repo / "output" / f"{args.infer_opt}_dim{embed_dim}_patch{patch_size}" / inst_name / "aln.tre"
        if not tre.exists():
            raise FileNotFoundError(f"Expected NeuralNJ output tree not found: {tre}")

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(tre, args.output)
        print(f"[neuralnj] wrote tree -> {args.output}")
    finally:
        shutil.rmtree(inst_dir, ignore_errors=True)
        os.unlink(cfg_path)


if __name__ == "__main__":
    main()
