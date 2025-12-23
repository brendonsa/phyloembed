rule simulate_species_tree:
    output:
        "results/species_trees/{n}/species_{n}taxa_{i}.nwk"
    params:
        height=config["tree_height_years"],
        lambda_=config["speciation_rate"],
        mu=config["extinction_rate"],
    conda:
        "phyloembed"
    shell:
        r"""
        mkdir -p results/species_trees/{wildcards.n}
        Rscript scripts/simulate_one_species_tree.R \
          --n_tips {wildcards.n} \
          --tree_height {params.height} \
          --lambda {params.lambda_} \
          --mu {params.mu} \
          --seed {wildcards.i} \
          --out {output}
        """

rule simphy_gene_tree:
    input:
        stree="results/species_trees/{n}/species_{n}taxa_{i}.nwk"
    output:
        done="results/simphy/{n}/dl={dl}/ps={ps}/rep={i}/done.txt"
    params:
        outdir=lambda wc: f"results/simphy/{wc.n}/dl={wc.dl}/ps={wc.ps}/rep={wc.i}",
        simphy=config["simphy_bin"],
        simphy_params = "config/simphy.params",
        seed=22,
    conda:
        "phyloembed"
    shell:
        r"""
        mkdir -p {params.outdir}

        TREE_STR="$(tr -d '\n' < {input.stree})"

        {params.simphy} \
          -S "$TREE_STR" \
          -I {params.simphy_params}
        echo "OK" > {output.done}
        """

