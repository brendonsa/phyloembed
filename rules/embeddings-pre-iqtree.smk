wildcard_constraints:
    dataset = "[^/]+",
    bootkind = "|".join(BOOTKINDS),
    mode = "|".join(MODES),
    model = "|".join(config["models"]),
    metric = "|".join(config["metrics"])

rule extract_embeddings:
    input:
        fasta = lambda wc: (
            f"results/{wc.dataset}/aligned.fasta"
            if wc.mode == "alignment"
            else f"results/{wc.dataset}/translated.fasta"
        )
    output:
        csv = "results/{dataset}/embeddings/{model}_{mode}.csv"
    conda:
        "phyloembed"
    resources:
        gpu_slots = 1
    shell:
        "bash scripts/extract_from_model.sh {wildcards.dataset} {wildcards.model} {input.fasta} {output.csv}"


rule compute_distance_matrix:
    input:
        csv = rules.extract_embeddings.output.csv
    output:
        dist = "results/{dataset}/distances_{mode}_{model}_{metric}.csv"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/compute_distance_matrix_phylip.py \
            --input {input.csv} \
            --output {output.dist} \
            --metric {wildcards.metric}
        """

rule build_nj_tree:
    input:
        dist = rules.compute_distance_matrix.output.dist
    output:
        tree = "results/{dataset}/tree_nj_{mode}_{model}_{metric}.nwk"
    conda:
        "phyloembed"
    shell:
        "fastme -i {input.dist} -o {output.tree} --nni --spr -q"


rule generate_bootstrap_inputs:
    input:
        aligned   = "results/{dataset}/aligned.fasta",
        unaligned = "results/{dataset}/translated.fasta"
    output:
        expand(
            "results/{{dataset}}/bootstrap/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        ),
        expand(
            "results/{{dataset}}/pseudo_bootstrap/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        )
    params:
        n = config.get("bootstrap_reps", 100)
    conda:
        "phyloembed"
    shell:
        """
        python scripts/generate_bootstrap_inputs.py \
            --aligned {input.aligned} \
            --unaligned {input.unaligned} \
            --output-base results/{wildcards.dataset} \
            --num {params.n}
        """

rule generate_pmsa_bootstrap_inputs:
    input:
        aligned = "results/{dataset}/aligned.fasta"
    output:
        expand(
            "results/{{dataset}}/pmsabootstrapw/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        ),
        expand(
            "results/{{dataset}}/pmsabootstrapwithout/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        )
    params:
        n    = config.get("bootstrap_reps", 100),
        frac = config.get("pmsa_scramble_frac", 0.05)
    conda:
        "phyloembed"
    shell:
        """
        python scripts/generate_pmsa_bootstrap_inputs.py \
          --aligned {input.aligned} \
          --output-base results/{wildcards.dataset} \
          --num {params.n} \
          --fraction {params.frac}
        """

rule generate_pmsa_bootstrap_inputs_frac50:
    input:
        aligned = "results/{dataset}/aligned.fasta"
    output:
        expand(
            "results/{{dataset}}/pmsabootstrapw_frac50/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        ),
        expand(
            "results/{{dataset}}/pmsabootstrapwithout_frac50/{rep}/boot.fasta",
            rep=[str(i) for i in range(1, config["bootstrap_reps"] + 1)],
        )
    params:
        n      = config.get("bootstrap_reps", 100),
        frac   = 0.5,
        suffix = "_frac50"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/generate_pmsa_bootstrap_inputs.py \
          --aligned {input.aligned} \
          --output-base results/{wildcards.dataset} \
          --num {params.n} \
          --fraction {params.frac} \
          --suffix {params.suffix}
        """


rule init_gpu_slots:
    output:
        "gpu_slots.txt"
    params:
        num_gpus    = config["num_gpus"],
        jobs_per_gpu = config["jobs_per_gpu"]
    shell:
        """
        rm -f gpu_slots.txt
        for gpu in $(seq 0 $(( {params.num_gpus} - 1 ))); do
            for slot in $(seq 1 {params.jobs_per_gpu}); do
                echo $gpu >> gpu_slots.txt
            done
        done
        echo "Initialized gpu_slots.txt"
        """


rule embed_bootstrap_rep:
    input:
        fasta = "results/{dataset}/{bootkind}/{rep}/boot.fasta",
    output:
        csv = "results/{dataset}/{bootkind}/{rep}/embeddings/{model}.csv"
    conda:
        "phyloembed"
    resources:
        gpu_slots = 1
    shell:
        "bash scripts/extract_from_model.sh {wildcards.dataset} {wildcards.model} {input.fasta} {output.csv}"


rule compute_bootstrap_distance_matrix:
    input:
        csv = rules.embed_bootstrap_rep.output.csv
    output:
        dist = "results/{dataset}/{bootkind}/{rep}/distances_{model}_{metric}.csv"
    conda:
        "phyloembed"
    shell:
        """
        python scripts/compute_distance_matrix_phylip.py \
            --input {input.csv} \
            --output {output.dist} \
            --metric {wildcards.metric}
        """


rule build_bootstrap_nj_embedding_tree:
    input:
        dist = rules.compute_bootstrap_distance_matrix.output.dist
    output:
        tree = "results/{dataset}/{bootkind}/{rep}/tree_nj_embedding_{model}_{metric}.nwk"
    conda:
        "phyloembed"
    shell:
        "fastme -i {input.dist} -o {output.tree} --nni --spr -q"


rule build_bootstrap_nj_tree:
    input:
        msa = "results/{dataset}/bootstrap/{rep}/boot.fasta"
    output:
        tree = "results/{dataset}/bootstrap/{rep}/tree_nj_msa.nwk"
    conda:
        "phyloembed"
    shell:
        "python scripts/build_nj_tree.py -i {input.msa} -o {output.tree} --from-msa"


rule build_base_bootstrap_tree:
    input:
        msa = "results/{dataset}/aligned.fasta"
    output:
        tree = "results/{dataset}/bootstrap/tree_nj.nwk"
    conda:
        "phyloembed"
    shell:
        "python scripts/build_nj_tree.py -i {input.msa} -o {output.tree} --from-msa"


rule build_bootstrap_with_support_trees:
    input:
        base_tree = rules.build_base_bootstrap_tree.output.tree,
        bootstrap_trees = expand(
            rules.build_bootstrap_nj_tree.output.tree,
            dataset="{dataset}",
            rep=[str(i) for i in range(1, config['bootstrap_reps'] + 1)],
        )
    output:
        supported = "results/{dataset}/tree_nj_msa_with_support.nwk",
    conda:
        "phyloembed"
    shell:
        """
        python scripts/bootstrap_nj_trees.py \
          --base-tree {input.base_tree} \
          --trees {input.bootstrap_trees} \
          --out-supported {output.supported} \
          --cutoff 0.5
        """


rule build_bootstrap_embedding_with_support_trees:
    input:
        base_tree = lambda wc: rules.build_nj_tree.output.tree.format(
            dataset=wc.dataset,
            mode="alignment",
            model=wc.model,
            metric=wc.metric,
        ),
        bootstrap_trees = expand(
            rules.build_bootstrap_nj_embedding_tree.output.tree,
            dataset="{dataset}",
            bootkind="{bootkind}",
            model="{model}",
            metric="{metric}",
            rep=[str(i) for i in range(1, config['bootstrap_reps'] + 1)],
        )
    output:
        supported = "results/{dataset}/{bootkind}_tree_nj_embedding_{model}_{metric}_with_support.nwk",
    conda:
        "phyloembed"
    shell:
        """
        python scripts/bootstrap_nj_trees.py \
          --base-tree {input.base_tree} \
          --trees {input.bootstrap_trees} \
          --out-supported {output.supported} \
          --cutoff 0.5
        """
