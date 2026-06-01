rule sanitize_fasta_labels:
    input:
        fasta = "data/{dataset}.fasta"
    output:
        clean = "results/{dataset}/clean.fasta",
        map   = "results/{dataset}/clean.labels.tsv"
    conda:
        "phyloembed"
    shell:
        "python scripts/sanitize_fasta.py -i {input.fasta} -o {output.clean} --map {output.map}"


rule align_nucleotides:
    input:
        fasta = "results/{dataset}/clean.fasta"
    output:
        aln = "results/{dataset}/aligned_nucleotides.fasta"
    conda:
        "phyloembed"
    shell:
        "mafft --auto {input.fasta} > {output.aln}"

rule translate_nucleotides:
    input:
        fasta = "results/{dataset}/clean.fasta"
    output:
        prot = "results/{dataset}/translated.fasta"
    conda:
        "phyloembed"
    shell:
        "python scripts/translate_fasta.py -i {input.fasta} -o {output.prot}"




rule align_proteins:
    input:
        prot = "results/{dataset}/translated.fasta"
    output:
        aln = "results/{dataset}/aligned.fasta"
    conda:
        "phyloembed"
    shell:
        "mafft --auto {input.prot} > {output.aln}"


rule infer_tree_raxml_nt:
    input:
        aln = "results/{dataset}/aligned_nucleotides.fasta"
    output:
        tree = "results/{dataset}/tree_raxml_nucleotides.nwk"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/raxml_nt
        raxml-ng --msa {input.aln} \
                 --model GTR+G \
                 --threads 2 \
                 --seed 42 \
                 --prefix results/{wildcards.dataset}/raxml_nt/{wildcards.dataset}_nt \
                 --redo
        cp results/{wildcards.dataset}/raxml_nt/{wildcards.dataset}_nt.raxml.bestTree {output.tree}
        """


rule infer_tree_raxml_aa:
    input:
        aln = "results/{dataset}/aligned.fasta"
    output:
        tree = "results/{dataset}/tree.nwk"
    conda: "phyloembed"
    threads: 2
    shell:
        """
        mkdir -p results/{wildcards.dataset}/raxml_aa
        raxml-ng --msa {input.aln} \
                 --model LG+G --threads {threads} --seed 42 \
                 --prefix results/{wildcards.dataset}/raxml_aa/{wildcards.dataset}_aa --redo
        cp results/{wildcards.dataset}/raxml_aa/{wildcards.dataset}_aa.raxml.bestTree {output.tree}
        """


rule use_external_alignments:
    input:
        ext = lambda wc: f"data/{wc.dataset}_aligned.fasta",
        map = "results/{dataset}/clean.labels.tsv"
    output:
        alignment = "results/{dataset}/aligned.fasta"
    conda:
        "phyloembed"
    shell:
        "python scripts/sanitize_fasta.py -i {input.ext} -o {output.alignment} --map {input.map}"


rule use_external_tree:
    input:
        ext = lambda wc: f"data/{wc.dataset}.nwk",
        map = "results/{dataset}/clean.labels.tsv"
    output:
        tree = "results/{dataset}/tree.nwk"
    conda:
        "phyloembed"
    shell:
        "python scripts/sanitize_tree_from_fasta.py --tree {input.ext} --map {input.map} -o {output.tree}"
