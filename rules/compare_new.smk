RESULTS_DIR = globals().get("RESULTS_DIR", "results")

rule compare_gwindows_tree:
    input:
        ref = RESULTS_DIR + "/{dataset}/tree.nwk",
        t2  = RESULTS_DIR + "/{dataset}/tree_nj_gwindows_{model}_{metric}_w{winsize}_ov{overlap}_s{sigma}.nwk"
    output:
        csv = RESULTS_DIR + "/{dataset}/distances_new/gwindows_{model}_{metric}_w{winsize}_ov{overlap}_s{sigma}.csv"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p $(dirname {output.csv})
        Rscript scripts/compare_trees.R --t1 {input.ref} --t2 {input.t2} --out {output.csv}
        """

rule compare_phyla_tree:
    input:
        ref = RESULTS_DIR + "/{dataset}/tree.nwk",
        t2  = RESULTS_DIR + "/{dataset}/phyla.nwk"
    output:
        csv = RESULTS_DIR + "/{dataset}/distances_new/phyla.csv"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p $(dirname {output.csv})
        Rscript scripts/compare_trees.R --t1 {input.ref} --t2 {input.t2} --out {output.csv}
        """

rule compare_neuralnj_tree:
    input:
        ref = RESULTS_DIR + "/{dataset}/tree.nwk",
        t2  = RESULTS_DIR + "/{dataset}/neuralnj.nwk"
    output:
        csv = RESULTS_DIR + "/{dataset}/distances_new/neuralnj.csv"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p $(dirname {output.csv})
        Rscript scripts/compare_trees.R --t1 {input.ref} --t2 {input.t2} --out {output.csv}
        """

def new_distance_inputs(wc):
    files = [
        f"{RESULTS_DIR}/{wc.dataset}/distances_new/gwindows_{model}_{metric}_w{winsize}_ov{overlap}_s{sigma}.csv"
        for model in config["models_concat"]
        for metric in config["metrics"]
        for winsize in config["gwindows"]["sizes"]
        for overlap in config["gwindows"]["overlaps"]
        for sigma in config["gwindows"]["sigmas"]
    ]
    if config.get("gwindows", {}).get("include_win1", True):
        files += [
            f"{RESULTS_DIR}/{wc.dataset}/distances_new/gwindows_{model}_{metric}_w1_ov0_s1.csv"
            for model in config["models_concat"]
            for metric in config["metrics"]
        ]
    if config.get("run_phyla", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/distances_new/phyla.csv")
    if config.get("run_neuralnj", True):
        files.append(f"{RESULTS_DIR}/{wc.dataset}/distances_new/neuralnj.csv")
    return files

rule aggregate_new_distances:
    input:
        new_distance_inputs
    output:
        csv = RESULTS_DIR + "/{dataset}/tree_distances_vs_ref_new.csv"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/aggregate_new_distances.py --inputs {input} --out {output.csv}
        """
