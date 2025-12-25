import os

sim_parse_cfg = config.get("sim_parse", {})
sim_parse_in = sim_parse_cfg.get("input_dir", "simulations/raw")
sim_parse_out = sim_parse_cfg.get(
    "output_dir",
    globals().get("RESULTS_DIR", "simulations/results"),
)

rule sim_parse_one:
    input:
        tree=os.path.join(sim_parse_in, "{dataset}", "gene_tree.nwk"),
        phy=os.path.join(sim_parse_in, "{dataset}", "alignment_TRUE.phy")
    output:
        tree=os.path.join(sim_parse_out, "{dataset}", "tree.nwk"),
        aligned=os.path.join(sim_parse_out, "{dataset}", "aligned.fasta")
    conda:
        "phyloembed"
    shell:
        r"""
        mkdir -p $(dirname {output.tree})
        cp {input.tree} {output.tree}
        python scripts/convert_phylip_to_fasta.py \
            --input {input.phy} \
            --output {output.aligned}
        """
