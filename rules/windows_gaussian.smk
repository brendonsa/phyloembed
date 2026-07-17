RESULTS_DIR = globals().get("RESULTS_DIR", "results")

wildcard_constraints:
    winsize = r"\d+",
    overlap = r"\d+",
    sigma   = r"\d+"

rule gaussian_window_distances:
    input:
        emb = RESULTS_DIR + "/{dataset}/embeddings/{model}_raw.npz"
    output:
        joined = RESULTS_DIR + "/{dataset}/gwindows/model={model}__metric={metric}__w={winsize}__ov={overlap}__s={sigma}/joined.phy",
        meta   = RESULTS_DIR + "/{dataset}/gwindows/model={model}__metric={metric}__w={winsize}__ov={overlap}__s={sigma}/windows.json"
    conda:
        "phyloembed"
    params:
        joined_only = int(config.get("gwindows", {}).get("joined_only", False))
    shell:
        r"""
        python scripts/build_gaussian_window_distances.py \
            --embeddings {input.emb} \
            --outdir $(dirname {output.joined}) \
            --window-size {wildcards.winsize} \
            --overlap {wildcards.overlap} \
            --sigma {wildcards.sigma} \
            --metric {wildcards.metric} \
            --reduce mean \
            $( [ "{params.joined_only}" = "1" ] && echo "--joined-only" || true )
        """

rule gaussian_window_phy_to_tree:
    input:
        dist = RESULTS_DIR + "/{dataset}/gwindows/model={model}__metric={metric}__w={winsize}__ov={overlap}__s={sigma}/joined.phy"
    output:
        tree = RESULTS_DIR + "/{dataset}/tree_nj_gwindows_{model}_{metric}_w{winsize}_ov{overlap}_s{sigma}.nwk"
    conda:
        "phyloembed"
    shell:
        r"""
        set -euo pipefail
        fastme -i {input.dist} -o {output.tree} --nni --spr -q  -T 8 >/dev/null
        """
