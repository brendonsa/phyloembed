MODELS_ALL = sorted(set(config["models"] + config.get("models_concat", [])))

wildcard_constraints:
    dataset = "[^/]+",
    bootkind = "|".join(BOOTKINDS),
    mode = "|".join(MODES),
    model = "|".join(MODELS_ALL),
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


rule extract_embeddings_raw:
    input:
        fasta = "results/{dataset}/aligned.fasta"
    output:
        npz = "results/{dataset}/embeddings/{model}_raw.npz"
    conda:
        "phyloembed"
    resources:
        gpu_slots = 1
    shell:
        "bash scripts/extract_from_model.sh {wildcards.dataset} {wildcards.model} {input.fasta} {output.npz} --pooling concat"


rule site_weights:
    input:
        msa="results/{dataset}/aligned.fasta"
    output:
        weights="results/{dataset}/site_weights_{wkind}.npy"
    params:
        wkind=lambda wc: wc.wkind
    conda:
        "phyloembed"
    shell:
        r"""
        python scripts/site_weights.py \
            --input {input.msa} \
            --output {output.weights} \
            --kind {params.wkind}
        """

def windows_param_folder(wc):
    return f"model={wc.model}__metric={wc.metric}__w={wc.winsize}__ov={wc.overlap}"

rule window_distances:
    input:
        emb = rules.extract_embeddings_raw.output.npz
    output:
        joined="results/{dataset}/windows/model={model}__metric={metric}__w={winsize}__ov={overlap}/joined.phy",
        meta="results/{dataset}/windows/model={model}__metric={metric}__w={winsize}__ov={overlap}/windows.json"
    conda:
        "phyloembed"
    shell:
        r"""
        python scripts/build_window_distances.py \
            --embeddings {input.emb} \
            --outdir $(dirname {output.joined}) \
            --window-size {wildcards.winsize} \
            --overlap {wildcards.overlap} \
            --metric {wildcards.metric} \
            --reduce mean
        """

rule sorted_coverage_distance:
    input:
        emb = rules.extract_embeddings_raw.output.npz,
        w = rules.site_weights.output.weights
    output:
        phy="results/{dataset}/sorted_cov/model={model}__metric={metric}__wkind={wkind}__inc0={inc0}__cov={coverage}pct__dir={direction}__norm={norm}/joined.phy",
        meta="results/{dataset}/sorted_cov/model={model}__metric={metric}__wkind={wkind}__inc0={inc0}__cov={coverage}pct__dir={direction}__norm={norm}/joined.json"
    conda:
        "phyloembed"
    shell:
        r"""
        set -euo pipefail
        python scripts/build_sorted_coverage_distance.py \
            --embeddings {input.emb} \
            --weights {input.w} \
            --out {output.phy} \
            --meta {output.meta} \
            --coverage $(python - <<'PY'
print(float("{wildcards.coverage}")/100.0)
PY
) \
            --metric {wildcards.metric} \
            --direction {wildcards.direction} \
            --normalize {wildcards.norm} \
            $( [ "{wildcards.inc0}" = "1" ] && echo "--include-zero" || true )
        """

rule window_phy_to_tree:
    input:
        dist="results/{dataset}/windows/model={model}__metric={metric}__w={winsize}__ov={overlap}/joined.phy"
    output:
        tree="results/{dataset}/tree_nj_windows_{model}_{metric}_w{winsize}_ov{overlap}.nwk"
    conda:
        "phyloembed"
    shell:
        r"""
        set -euo pipefail
        fastme -i {input.dist} -o {output.tree} --nni --spr -q  -T 1 >/dev/null
        """

rule sorted_cov_phy_to_tree:
    input:
        dist="results/{dataset}/sorted_cov/model={model}__metric={metric}__wkind={wkind}__inc0={inc0}__cov={coverage}pct__dir={direction}__norm={norm}/joined.phy"
    output:
        tree="results/{dataset}/tree_nj_sortedcov_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{direction}_{norm}.nwk"
    conda:
        "phyloembed"
    shell:
        r"""
        set -euo pipefail
        fastme -i {input.dist} -o {output.tree} --nni --spr -q  -T 1 >/dev/null
        """

def subset_param_folder(wc):
    return (
        f"model={wc.model}__metric={wc.metric}"
        f"__cov={wc.coverage}__n={wc.nsubsets}"
        f"__seed={wc.seed}__wkind={wc.wkind}"
    )

rule weighted_subset_consensus:
    input:
        emb = rules.extract_embeddings_raw.output.npz,
        w = rules.site_weights.output.weights
    output:
        majority="results/{dataset}/consensus_weighted_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{nsubsets}_seed{seed}_majority.nwk",
        greedy="results/{dataset}/consensus_weighted_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{nsubsets}_seed{seed}_greedy.nwk",
        meta="results/{dataset}/subsets_consensus/meta_{model}_{metric}_{wkind}_inc0{inc0}_{coverage}pct_{nsubsets}_seed{seed}.json"
    conda:
        "phyloembed"
    shell:
        r"""
        set -euo pipefail

        WORKDIR="results/{wildcards.dataset}/subsets_consensus/work_model={wildcards.model}__metric={wildcards.metric}__wkind={wildcards.wkind}__inc0={wildcards.inc0}__cov={wildcards.coverage}pct__n={wildcards.nsubsets}__seed={wildcards.seed}"
        mkdir -p "$WORKDIR"

        python scripts/build_weighted_subset_phylip.py \
            --embeddings {input.emb} \
            --weights {input.w} \
            --outdir "$WORKDIR" \
            --n-subsets {wildcards.nsubsets} \
            --coverage $(python - <<'PY'
cov = float("{wildcards.coverage}")/100.0
print(cov)
PY
) \
            --metric {wildcards.metric} \
            --seed {wildcards.seed} \
            $( [ "{wildcards.inc0}" = "1" ] && echo "--include-zero" || true )

        cp "$WORKDIR/subsets.json" {output.meta}

        rm -f "$WORKDIR"/distance_*.nwk "$WORKDIR/subset_trees.nwk" || true
        for phy in "$WORKDIR"/distance_*.phy; do
            base="$(basename "$phy" .phy)"
            nwk="$WORKDIR/${{base}}.nwk"
            fastme -i "$phy" -o "$nwk" --nni --spr   -T 1 >/dev/null
            cat "$nwk" >> "$WORKDIR/subset_trees.nwk"
            echo >> "$WORKDIR/subset_trees.nwk"
        done

        iqtree2 -con "$WORKDIR/subset_trees.nwk" -minsup 0.5 -pre "$WORKDIR/majority" -T 1 >/dev/null 2>&1
        cp "$WORKDIR/majority.contree" {output.majority}

        iqtree2 -con "$WORKDIR/subset_trees.nwk" -minsup 0 -pre "$WORKDIR/greedy" -T 1 >/dev/null 2>&1
        cp "$WORKDIR/greedy.contree" {output.greedy}
        """

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
        "fastme -i {input.dist} -o {output.tree} --nni --spr -q  -T 1 >/dev/null" 


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
        scripts/iqtree_support.sh \
          {input.base_tree} \
          {input.bootstrap_trees} \
          {output.supported}
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
        scripts/iqtree_support.sh \
          {input.base_tree} \
          {input.bootstrap_trees} \
          {output.supported}
        """
