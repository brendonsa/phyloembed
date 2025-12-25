RESULTS_DIR = globals().get("RESULTS_DIR", "results")

rule phyloformer_prepare_aln:
    input:
        msa = RESULTS_DIR + "/{dataset}/aligned.fasta"
    output:
        aln = RESULTS_DIR + "/{dataset}/phyloformer/alns/aln.fasta"
    params:
        aln_dir = RESULTS_DIR + "/{dataset}/phyloformer/alns"
    shell:
        """
        rm -rf {params.aln_dir}
        mkdir -p {params.aln_dir}
        cp {input.msa} {output.aln}
        """

rule phyloformer_infer_matrix:
    input:
        aln = rules.phyloformer_prepare_aln.output.aln
    output:
        mat = RESULTS_DIR + "/{dataset}/phyloformer/pf_matrices/aln.phy"
    params:
        model = config.get("phyloformer_model", "Phyloformer/models/pf.ckpt"),
        mats_dir = RESULTS_DIR + "/{dataset}/phyloformer/pf_matrices"
    conda:
        "phyloformer"
    shell:
        r"""
        mkdir -p {params.mats_dir}
        WORKDIR="$(mktemp -d)"
        cp {input.aln} "$WORKDIR/aln.fasta"
        python Phyloformer/infer_alns.py -o {params.mats_dir} {params.model} "$WORKDIR"
        rm -rf "$WORKDIR"
        """

rule phyloformer_fastme_tree:
    input:
        mat = rules.phyloformer_infer_matrix.output.mat
    output:
        tree = RESULTS_DIR + "/{dataset}/phyloformer.nwk"
    params:
    conda:
        "phyloformer" 
    shell:
        r"""
        fastme -i {input.mat} -o {output.tree} --nni --spr
        """

rule phyloformer:
    input:
        rules.phyloformer_fastme_tree.output.tree
