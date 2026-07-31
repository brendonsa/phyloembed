RESULTS_DIR = globals().get("RESULTS_DIR", "results")

def new_tree_paths(wc):
    files = [
        f"{RESULTS_DIR}/{wc.dataset}/tree_nj_gwindows_{model}_{metric}_w{winsize}_ov{overlap}_s{sigma}.nwk"
        for model in config["models_concat"]
        for metric in config["metrics"]
        for winsize in config["gwindows"]["sizes"]
        for overlap in config["gwindows"]["overlaps"]
        for sigma in config["gwindows"]["sigmas"]
    ]
    if config.get("gwindows", {}).get("include_win1", True):
        files += [
            f"{RESULTS_DIR}/{wc.dataset}/tree_nj_gwindows_{model}_{metric}_w1_ov0_s1.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
        ]
    if config.get("run_phyla", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/phyla.nwk")
    if config.get("run_neuralnj", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/neuralnj.nwk")
    if config.get("run_iqtree", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/tree_iqtree.nwk")
    if config.get("run_fasttree", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/tree_fasttree.nwk")
    return files


rule build_iqtree_tree:
    input:
        aln = RESULTS_DIR + "/{dataset}/aligned.fasta"
    output:
        tree = RESULTS_DIR + "/{dataset}/tree_iqtree.nwk"
    conda:
        "phyloembed"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/{wildcards.dataset}/iqtree
        iqtree2 -s {input.aln} -m LG+G -T {threads} \
            --prefix results/{wildcards.dataset}/iqtree/{wildcards.dataset} -redo >/dev/null
        cp results/{wildcards.dataset}/iqtree/{wildcards.dataset}.treefile {output.tree}
        """


rule build_fasttree_tree:
    input:
        aln = RESULTS_DIR + "/{dataset}/aligned.fasta"
    output:
        tree = RESULTS_DIR + "/{dataset}/tree_fasttree.nwk"
    conda:
        "phyloembed"
    shell:
        "FastTree -lg -quiet {input.aln} > {output.tree}"

rule tree_distances_vs_ref_new:
    input:
        old_csv = RESULTS_DIR + "/{dataset}/tree_distances_vs_ref.csv",
        new_trees = new_tree_paths
    output:
        csv = RESULTS_DIR + "/{dataset}/tree_distances_vs_ref_new.csv"
    conda:
        "phyloembed"
    params:
        dataset_dir = RESULTS_DIR + "/{dataset}"
    shell:
        """
        Rscript scripts/compare_all_trees.r \
            --dataset-dir {params.dataset_dir} \
            --out {output.csv}
        """
