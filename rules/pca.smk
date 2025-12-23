# PCA rules

rule extract_embeddings_pca:
    input:
        fasta = lambda wc: (
            f"results/{wc.dataset}/aligned.fasta"
            if wc.mode == "alignment"
            else f"results/{wc.dataset}/translated.fasta"
        )
    output:
        csv = "results/{dataset}/embeddings/{model}_{mode}_pca.csv"
    conda:
        "phyloembed"
    resources:
        gpu_slots = 1
    shell:
        """
        bash scripts/extract_from_model.sh \
            {wildcards.dataset} \
            {wildcards.model} \
            {input.fasta} \
            {output.csv} \
            --pooling pca
        """

rule compute_distance_matrix_pca:
    input:
        csv = "results/{dataset}/embeddings/{model}_{mode}_pca.csv"
    output:
        dist = "results/{dataset}/distances_{mode}_{model}_pca_{metric}.csv"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/compute_distance_matrix.py \
            --input {input.csv} \
            --output {output.dist} \
            --metric {wildcards.metric}
        """

rule build_nj_tree_pca:
    input:
        dist = "results/{dataset}/distances_{mode}_{model}_pca_{metric}.csv"
    output:
        tree = "results/{dataset}/tree_nj_{mode}_{model}_pca_{metric}.nwk"
    conda:
        "phyloembed"
    shell:
        "python scripts/build_nj_tree.py -i {input.dist} -o {output.tree}"

rule embed_bootstrap_rep_pca:
    input:
        fasta = "results/{dataset}/{bootkind}/{rep}/boot.fasta",
        slots = ancient("gpu_slots.txt")
    output:
        csv = "results/{dataset}/{bootkind}/{rep}/embeddings/{model}_pca.csv"
    conda:
        "phyloembed"
    resources:
        gpu_slots=1
    shell:
        """
        bash scripts/extract_from_model.sh \
            {wildcards.dataset} \
            {wildcards.model} \
            {input.fasta} \
            {output.csv} \
            --pooling pca
        """

rule compute_bootstrap_distance_matrix_pca:
    input:
        csv = "results/{dataset}/{bootkind}/{rep}/embeddings/{model}_pca.csv"
    output:
        dist = "results/{dataset}/{bootkind}/{rep}/distances_{model}_pca_{metric}.csv"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/compute_distance_matrix.py \
            --input {input.csv} \
            --output {output.dist} \
            --metric {wildcards.metric}
        """

rule build_bootstrap_nj_embedding_tree_pca:
    input:
        dist = "results/{dataset}/{bootkind}/{rep}/distances_{model}_pca_{metric}.csv"
    output:
        tree = "results/{dataset}/{bootkind}/{rep}/tree_nj_embedding_{model}_pca_{metric}.nwk"
    conda:
        "phyloembed"
    shell:
        "python scripts/build_nj_tree.py -i {input.dist} -o {output.tree}"

rule build_bootstrap_embedding_consensus_trees_pca:
    input:
        base_tree = "results/{dataset}/tree_nj_alignment_{model}_pca_{metric}.nwk",
        bootstrap_trees = expand(
            "results/{{dataset}}/{{bootkind}}/{rep}/tree_nj_embedding_{{model}}_pca_{{metric}}.nwk",
            rep=[str(i) for i in range(1, config['bootstrap_reps'] + 1)],
        )
    output:
        supported = "results/{dataset}/{bootkind}_tree_nj_embedding_{model}_pca_{metric}_with_supports.nwk",
        consensus = "results/{dataset}/{bootkind}_tree_nj_embedding_{model}_pca_{metric}_consensus.nwk"
    conda: "phyloembed"
    shell:
        """
        python scripts/bootstrap_nj_trees.py \
          --base-tree {input.base_tree} \
          --trees {input.bootstrap_trees} \
          --out-supported {output.supported} \
          --out-consensus {output.consensus} \
          --cutoff 0.5
        """
