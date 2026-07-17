RESULTS_DIR = globals().get("RESULTS_DIR", "results")

# NeuralNJ only models nucleotides (DNA_WITH_GAP). Feed it the nucleotide alignment,
# not the protein MSA. Repo assumed checked out as a sibling dir (config: neuralnj_repo).

rule neuralnj_prepare_phy:
    input:
        aln = RESULTS_DIR + "/{dataset}/aligned_nucleotides.fasta"
    output:
        phy = RESULTS_DIR + "/{dataset}/neuralnj/aln.phy"
    conda:
        "phyloembed"
    shell:
        """
        mkdir -p $(dirname {output.phy})
        python scripts/fasta_to_phylip.py -i {input.aln} -o {output.phy}
        """

rule neuralnj_tree:
    input:
        phy = rules.neuralnj_prepare_phy.output.phy
    output:
        tree = RESULTS_DIR + "/{dataset}/neuralnj.nwk"
    params:
        repo       = config.get("neuralnj_repo", "NeuralNJ"),
        infer_opt  = config.get("neuralnj_infer_opt", "Argmax"),
        checkpoint = config.get("neuralnj_checkpoint", "")
    conda:
        "neuralnj"
    resources:
        gpu_slots = 1
    shell:
        r"""
        python scripts/run_neuralnj.py \
            --phy {input.phy} \
            --output {output.tree} \
            --repo {params.repo} \
            --infer-opt {params.infer_opt} \
            $( [ -n "{params.checkpoint}" ] && echo "--checkpoint {params.checkpoint}" || true )
        """

rule neuralnj:
    input:
        rules.neuralnj_tree.output.tree
