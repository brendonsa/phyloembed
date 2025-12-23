import os
import tempfile
import subprocess

CONDA = "/home/compute2/miniconda3/bin/conda"

def consensus_from_bootstrap_trees(tree_paths, out_tree, minsup=0.5):
    """
    tree_paths: list of .nwk paths (your 'trees' list)
    out_tree:   path to write consensus Newick
    minsup:     0.5 = majority rule (50%), 0 = extended MR
    """
    # 1) combine all bootstrap trees into one file
    with tempfile.NamedTemporaryFile(suffix=".tre", delete=False) as tmp:
        combined = tmp.name

    with open(combined, "w") as out:
        for p in tree_paths:
            with open(p) as f:
                txt = f.read()
                out.write(txt)
                if not txt.endswith("\n"):
                    out.write("\n")

    # 2) run iqtree2 consensus
    prefix = out_tree + ".tmp"
    subprocess.run(
        [
            CONDA, "run", "-n", "phyloembed",
            "iqtree2",
            "-t", combined,
            "-con",
            "-minsup", str(minsup),
            "-pre", prefix,
            "-nt", "1",
        ],
        check=True,
    )

    os.replace(prefix + ".contree", out_tree)

    # 3) clean up
    for ext in (".log", ".iqtree", ".trees"):
        f = prefix + ext
        if os.path.exists(f):
            os.remove(f)
    os.remove(combined)

    return out_tree