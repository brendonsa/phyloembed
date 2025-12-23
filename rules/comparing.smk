rule compare_trees:
    input:
        t1 = "results/{dataset}/tree_raxml_{type}.nwk",
        t2 = "results/{dataset}/tree_nj_{model}_{metric}.nwk"
    output:
        "results/{dataset}/distances/{type}_{model}_{metric}.csv"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/distances
        Rscript scripts/compare_trees.R \
          --t1 {input.t1} \
          --t2 {input.t2} \
          --out {output}
        """

rule compare_aa_nucleotides:
    input:
        t1 = "results/{dataset}/tree_raxml_nucleotides.nwk",
        t2 = "results/{dataset}/tree_raxml_aminoacids.nwk",
    output:
        "results/{dataset}/distances/nucleotides_aminoacids.csv"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p results/{wildcards.dataset}/distances
        Rscript scripts/compare_trees.R \
          --t1 {input.t1} \
          --t2 {input.t2} \
          --out {output}
        """