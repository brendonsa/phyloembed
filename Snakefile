configfile: "config/config.yaml"

MODES = ["alignment", "non_aligned"]
BOOTKINDS = [
#    "bootstrap",
#    "pseudo_bootstrap",
#    "pmsabootstrapw",
#    "pmsabootstrapwithout",
#    "pmsabootstrapw_frac50",
#    "pmsabootstrapwithout_frac50"
]
POOLINGS = [""]  # Empty string = no suffix = mean pooling (default)

include: "rules/sequence_analysis.smk"
include: "rules/embeddings.smk"
include: "rules/comparing.smk"
include: "rules/phyloformer.smk"

ruleorder: use_external_tree > infer_tree_raxml_aa
ruleorder: use_external_alignments > align_proteins

def subset_param_folder(model, metric, cov, nsubsets, seed, wkind):
    return (
        f"model={model}__metric={metric}"
        f"__cov={cov}__n={nsubsets}"
        f"__seed={seed}__wkind={wkind}"
    )

def window_param_folder(model, metric, winsize, overlap):
    return f"model={model}__metric={metric}__w={winsize}__ov={overlap}"



rule all:
    input:
        distances = expand(
            "results/{dataset}/tree_distances_vs_ref.csv",
            dataset=config["datasets"],
        ),
        raw_embeddings = expand(
            "results/{dataset}/embeddings/{model}_raw.npz",
            dataset=config["datasets"],
            model=config["models_concat"],
        )

        

rule tree_distances_vs_ref:
    input:
        ref="results/{dataset}/tree.nwk",
        msa="results/{dataset}/tree_nj_msa_with_support.nwk",
        plm=lambda wc: [
            f"results/{wc.dataset}/tree_nj_{mode}_{model}_{metric}.nwk"
            for mode in MODES
            for model in config["models"]
            for metric in config["metrics"]
        ]+ [
            f"results/{wc.dataset}/tree_nj_alignment_msapairformer_{metric}.nwk"
            for metric in config["metrics"]
        ],
        emb=lambda wc: [
            f"results/{wc.dataset}/{bootkind}_tree_nj_embedding_{model}_{metric}_with_support.nwk"
            for bootkind in BOOTKINDS
            for model in config["models"]
            for metric in config["metrics"]
        ],
        phyloformer="results/{dataset}/phyloformer.nwk",

        weighted_subset_majority=lambda wc: [
            f"results/{wc.dataset}/consensus_weighted_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{nsubsets}_seed{seed}_majority.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
            for wkind in config["weights"]
            for inc0 in config.get("include_zero", [0])
            for coverage in config["subsets"]["coverage_percentages"]
            for nsubsets in config["subsets"]["repetition_counts"]
            for seed in [config["subsets"]["seed"]]
        ],
        weighted_subset_greedy=lambda wc: [
            f"results/{wc.dataset}/consensus_weighted_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{nsubsets}_seed{seed}_greedy.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
            for wkind in config["weights"]
            for inc0 in config.get("include_zero", [0])
            for coverage in config["subsets"]["coverage_percentages"]
            for nsubsets in config["subsets"]["repetition_counts"]
            for seed in [config["subsets"]["seed"]]
        ],

        window_trees=lambda wc: [
            f"results/{wc.dataset}/tree_nj_windows_{model}_{metric}_w{winsize}_ov{overlap}.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
            for winsize in config["windows"]["sizes"]
            for overlap in config["windows"]["overlaps"]
        ] + [
            # explicit include of size1/ov0 even if not in lists
            f"results/{wc.dataset}/tree_nj_windows_{model}_{metric}_w1_ov0.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
        ] + [
            # explicit include of size0/ov0 even if not in lists
            f"results/{wc.dataset}/tree_nj_windows_{model}_{metric}_w0_ov0.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
        ],

        sorted_cov_trees=lambda wc: [
            f"results/{wc.dataset}/tree_nj_sortedcov_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{direction}_{norm}.nwk"
            for model in config["models_concat"]
            for metric in config["metrics"]
            for wkind in config["weights"]
            for inc0 in config["sorted_coverage"]["include_zero"]
            for coverage in config["sorted_coverage"]["coverage_percentages"]
            for direction in config["sorted_coverage"]["directions"]
            for norm in config["sorted_coverage"]["normalize"]
        ]
       
    conda:
        "phyloembed"
    params:
        dataset_dir = RESULTS_DIR + "/{dataset}"
    output:
        RESULTS_DIR + "/{dataset}/tree_distances_vs_ref.csv"
    shell:
        """
        Rscript scripts/compare_all_trees.r \
            --dataset-dir {params.dataset_dir} \
            --out {output}
        """