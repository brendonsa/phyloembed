rule simulate_tree:
    output:
        nwk=f"{trees_out}/L{{leaf}}_{{idx}}.nwk",
        nexus=f"{trees_out}/L{{leaf}}_{{idx}}.trees",
        times=f"{trees_out}/L{{leaf}}_{{idx}}.node_times.tsv"
    conda:
        "phyloembed"
    params:
        height=trees_cfg["height"],
        speciation_rate=trees_cfg["speciation_rate"],
        extinction_rate=trees_cfg["extinction_rate"],
        seed=lambda wc: int(trees_cfg["seed"]) + int(wc.leaf) * 1000 + int(wc.idx),
        script= "scripts/sim_species_tree.R"
    shell:
        (
            "Rscript {params.script} "
            "--leaf {wildcards.leaf} "
            "--height {params.height} "
            "--speciation_rate {params.speciation_rate} "
            "--extinction_rate {params.extinction_rate} "
            "--seed {params.seed} "
            "--out_nwk {output.nwk} "
            "--out_nexus {output.nexus} "
            "--out_times {output.times}"
        )


rule simphy_gene_trees:
    input:
        tree=f"{trees_out}/L{{leaf}}_{{idx}}.nwk"
    output:
        done=f"{sim_out}/L{{leaf}}_{{idx}}/ps{{ps}}_dl{{dl}}/done.txt"
    conda:
        "phyloembed"
    params:
        outdir=lambda wc: f"{sim_out}/L{wc.leaf}_{wc.idx}/ps{wc.ps}_dl{wc.dl}",
        binary=simphy_cfg["binary"],
        rs=simphy_cfg["rs"],
        parameters_file=simphy_cfg["parameters_file"]
    shell:
        (
            "mkdir -p {params.outdir} && "
            "TREE_STR=$(tr -d '\\n' < {input.tree}) && "
            "TREE_STR=\"${{TREE_STR%;}}\" && "
            "TREE_STR=$(printf \"%s\" \"$TREE_STR\" | sed -E 's/:[0-9.eE+-]+$//') && "
            "TREE_STR=\"$TREE_STR;\" && "
            "{params.binary} "
            "-S \"$TREE_STR\" "
            "-rs {params.rs} "
            "-I {params.parameters_file} "
            "-sp F:{wildcards.ps} "
            "-lb F:{wildcards.dl} "
            "-o {params.outdir} && "
            "touch {output.done}"
        )

