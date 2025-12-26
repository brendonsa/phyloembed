import os

sim_parse_cfg = config.get("sim_parse", {})
sim_parse_in = sim_parse_cfg.get("input_dir", "simulations/raw")
sim_parse_out = sim_parse_cfg.get(
    "output_dir",
    globals().get("RESULTS_DIR", "simulations/results"),
)

def sim_parse_seq_input(wc):
    base = os.path.join(sim_parse_in, wc.dataset)
    unaligned = os.path.join(base, "alignment.unaligned.fa")
    phy = os.path.join(base, "alignment_TRUE.phy")
    return unaligned if os.path.exists(unaligned) else phy

rule sim_parse_one:
    input:
        tree=os.path.join(sim_parse_in, "{dataset}", "gene_tree.nwk"),
        seq=sim_parse_seq_input
    output:
        tree=os.path.join(sim_parse_out, "{dataset}", "tree.nwk"),
        aligned=os.path.join(sim_parse_out, "{dataset}", "aligned.fasta"),
        unaligned=os.path.join(sim_parse_out, "{dataset}", "translated.fasta")
    conda:
        "phyloembed"
    shell:
        r"""
        mkdir -p $(dirname {output.tree})
        cp {input.tree} {output.tree}
        tmp_fasta="$(mktemp)"
        if [[ "{input.seq}" == *.phy ]]; then
            python scripts/convert_phylip_to_fasta.py \
                --input {input.seq} \
                --output "$tmp_fasta"
            cp "$tmp_fasta" {output.aligned}
            cp "$tmp_fasta" {output.unaligned}
        else
            cp {input.seq} {output.unaligned}
            mafft --auto {input.seq} > {output.aligned}
        fi
        rm -f "$tmp_fasta"
        """
