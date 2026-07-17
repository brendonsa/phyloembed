RESULTS_DIR = globals().get("RESULTS_DIR", "results")

rule phyla_prepare_aln:
    input:
        msa = RESULTS_DIR + "/{dataset}/aligned.fasta"
    output:
        aln = RESULTS_DIR + "/{dataset}/phyla/aln.fasta"
    shell:
        """
        mkdir -p $(dirname {output.aln})
        cp {input.msa} {output.aln}
        """

rule phyla_tree:
    input:
        aln = rules.phyla_prepare_aln.output.aln
    output:
        tree = RESULTS_DIR + "/{dataset}/phyla.nwk"
    params:
        model  = config.get("phyla_model", "phyla-beta"),
        device = config.get("phyla_device", "cuda:0")
    conda:
        "phyla"
    resources:
        gpu_slots = 1
    shell:
        r"""
        python scripts/run_phyla.py \
            --input {input.aln} \
            --output {output.tree} \
            --model {params.model} \
            --device {params.device}
        """

rule phyla:
    input:
        rules.phyla_tree.output.tree
