rule phyloformer_prepare_aln:
    input:
        msa = "results/{dataset}/aligned.fasta"
    output:
        aln = "results/{dataset}/phyloformer/alns/{dataset}.fasta"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/phyloformer/alns
        cp {input.msa} {output.aln}
        """

rule phyloformer_infer_matrix:
    input:
        aln = rules.phyloformer_prepare_aln.output.aln
    output:
        mat = "results/{dataset}/phyloformer/pf_matrices/{dataset}.phy"
    params:
        model = config.get("phyloformer_model", "Phyloformer/models/pf.ckpt"),
        mats_dir = "results/{dataset}/phyloformer/pf_matrices"
    conda:
        "phyloformer"
    shell:
        r"""
        mkdir -p {params.mats_dir}
        python Phyloformer/infer_alns.py -o {params.mats_dir} {params.model} $(dirname {input.aln})
        """

rule phyloformer_fastme_tree:
    input:
        mat = rules.phyloformer_infer_matrix.output.mat
    output:
        tree = "results/{dataset}/phyloformer.nwk"
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
