#!/bin/bash
# 
# Simulation Pipeline for MSAs and Gene Trees
# Based on Molloy and Warnow (2020) parameters
#
# This pipeline simulates:
# 1. Species trees
# 2. Gene trees (with duplication/loss and ILS)
# 3. DNA sequences, MSAs, and ML gene trees
# 4. Protein sequences, MSAs, and ML gene trees
#

set -e  # Exit on error

# ============================================================================
# USAGE FUNCTION
# ============================================================================

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Simulation pipeline for MSAs and gene trees (DNA and protein sequences).

OPTIONS:
    -n NUM        Number of replicates per configuration (default: 100)
    -l LEAVES     Comma-separated list of leaf counts (default: 12,20)
                  Example: -l 12,20,50
    -t TYPE       Sequence type to simulate: dna, protein, or both (default: both)
    -o DIR        Output directory (default: simulation_output)
    -s SEED       Random seed (default: 22)
    -r MAX        Retry if duplicate sequences found, up to MAX attempts (default: 0 = no retry)
    -f NUM        Filter gene trees: require exactly NUM leaves (default: 0 = no filtering)
    -m            Skip ML tree estimation (PhyML) - only generate alignments
    -i            Infer-only mode: only run PhyML on existing alignments (skip simulation)
    --dna-length NUM      DNA alignment length in base pairs (default: 1000)
    --protein-length NUM  Protein alignment length in amino acids (default: 333)
    --indel-model Indel preset: noindels, realistic, conservative, or highrate
                  noindels:     Disable indels (substitutions only)
                  realistic:    --indel 0.03,0.09 --indel-size "POW{1.7/50}" (default)
                  conservative: --indel 0.02,0.04 --indel-size "POW{1.7/20}"
                  highrate:     --indel 0.06,0.18 --indel-size "POW{1.7/50}"
    --indel       Enable indel simulation with rates INS,DEL (e.g., --indel 0.03,0.09)
    --indel-size  Indel length distribution: POW a max, GEO p, NB r p, or LAV a max
                  Example: --indel-size "POW 1.7 50" (default: POW 1.7 50)
    -h            Show this help message

EXAMPLES:
    # Run with default parameters (100 replicates, 12 and 20 leaves, both types)
    $0

    # Run 50 replicates with only 12 leaves, DNA only
    $0 -n 50 -l 12 -t dna

    # Run with 12, 20, and 50 leaves, protein only
    $0 -l 12,20,50 -t protein

    # Quick test: 10 replicates, 12 leaves, both types
    $0 -n 10 -l 12 -t both

    # Generate alignments only, skip ML tree estimation (faster)
    $0 -n 10 -l 12 -m

    # Infer ML trees from existing alignments (after running with -m)
    $0 -o test_output -i

    # No indels (substitutions only)
    $0 -n 10 -l 12 --indel-model noindels

    # Simulate with indels using preset model (recommended)
    $0 -n 10 -l 12 --indel-model realistic

    # Conservative indel model
    $0 -n 10 -l 12 --indel-model conservative

    # High indel rate for stress testing
    $0 -n 10 -l 12 --indel-model highrate

    # Custom indel parameters (manual specification)
    $0 -n 10 -l 12 --indel 0.03,0.09 --indel-size "POW 1.7 20"

    # Custom alignment lengths
    $0 -n 10 -l 12 --dna-length 500 --protein-length 200

    # Longer protein alignments with indels
    $0 -n 10 -l 12 --protein-length 500 --indel-model realistic

OUTPUT:
    The pipeline generates species trees, gene trees, and sequence alignments.
    By default, ML trees are estimated with PhyML. Use -m to skip ML estimation.
    Use -i to infer ML trees from existing alignments without re-running simulation.

EOF
    exit 1
}

# ============================================================================
# DEFAULT CONFIGURATION
# ============================================================================

# Number of replicates
NUM_REPLICATES=100

# Species tree parameters
TREE_HEIGHT=1800000337.5   # years
SPECIATION_RATE=1.8e-9     # events/year
EXTINCTION_RATE=0          # events/year
NUM_LEAVES=(12 20)         # Different leaf set sizes

# Gene tree parameters (SimPhy)
DUPLICATION_LOSS_RATES=(1e-10 2e-10 5e-10)
POPULATION_SIZES=(1e7 5e7)  # Low and medium ILS

# Sequence simulation parameters
ALIGNMENT_LENGTH_DNA=1000    # bp for DNA
ALIGNMENT_LENGTH_PROTEIN=333 # aa for proteins (~1000bp / 3)

# Random seed
RANDOM_SEED=22

# Output directories
OUTPUT_DIR="simulation_output"

# Run control
RUN_DNA=1
RUN_PROTEIN=1
RUN_ML_ESTIMATION=1  # 1 = run PhyML, 0 = skip ML tree estimation
INFER_ONLY=0         # 1 = only infer ML trees from existing alignments, skip simulation

# Duplicate handling
MAX_RETRIES=0  # 0 = no retry, >0 = retry up to MAX_RETRIES times

# Gene tree filtering
MIN_GENE_TREE_LEAVES=0  # 0 = no filtering, >0 = require exactly this many leaves
MAX_GENE_TREE_RETRIES=1000  # Maximum attempts to generate a tree with the required number of leaves

# Indel simulation parameters
ENABLE_INDELS=0
INDEL_RATES=""           # Format: "INS,DEL" e.g., "0.03,0.09"
INDEL_SIZE_DIST="POW{1.7/50}"  # Default: Zipfian with exponent 1.7, max length 50 (AliSim format)

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -n)
            NUM_REPLICATES=$2
            shift 2
            ;;
        -l)
            IFS=',' read -ra NUM_LEAVES <<< "$2"
            shift 2
            ;;
        -t)
            case ${2} in
                dna|DNA)
                    RUN_DNA=1
                    RUN_PROTEIN=0
                    ;;
                protein|PROTEIN)
                    RUN_DNA=0
                    RUN_PROTEIN=1
                    ;;
                both|BOTH)
                    RUN_DNA=1
                    RUN_PROTEIN=1
                    ;;
                *)
                    echo "Error: Invalid sequence type '${2}'. Use: dna, protein, or both"
                    usage
                    ;;
            esac
            shift 2
            ;;
        -o)
            OUTPUT_DIR=$2
            shift 2
            ;;
        -s)
            RANDOM_SEED=$2
            shift 2
            ;;
        -r)
            MAX_RETRIES=$2
            shift 2
            ;;
        -f)
            MIN_GENE_TREE_LEAVES=$2
            shift 2
            ;;
        -m)
            RUN_ML_ESTIMATION=0
            shift
            ;;
        -i)
            INFER_ONLY=1
            RUN_ML_ESTIMATION=1
            shift
            ;;
        --dna-length)
            ALIGNMENT_LENGTH_DNA=$2
            shift 2
            ;;
        --protein-length)
            ALIGNMENT_LENGTH_PROTEIN=$2
            shift 2
            ;;
        --indel-model)
            case ${2} in
                noindels)
                    ENABLE_INDELS=0
                    ;;
                realistic)
                    ENABLE_INDELS=1
                    INDEL_RATES="0.03,0.09"
                    INDEL_SIZE_DIST="POW{1.7/50}"
                    ;;
                conservative)
                    ENABLE_INDELS=1
                    INDEL_RATES="0.02,0.04"
                    INDEL_SIZE_DIST="POW{1.7/20}"
                    ;;
                highrate)
                    ENABLE_INDELS=1
                    INDEL_RATES="0.06,0.18"
                    INDEL_SIZE_DIST="POW{1.7/50}"
                    ;;
                *)
                    echo "Error: Invalid indel model '${2}'. Use: noindels, realistic, conservative, or highrate"
                    usage
                    ;;
            esac
            shift 2
            ;;
        --indel)
            ENABLE_INDELS=1
            INDEL_RATES=$2
            shift 2
            ;;
        --indel-size)
            INDEL_SIZE_DIST=$2
            shift 2
            ;;
        -h)
            usage
            ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            ;;
    esac
done

# Set directory paths based on output directory
SPECIES_TREES_DIR="${OUTPUT_DIR}/species_trees"
GENE_TREES_DIR="${OUTPUT_DIR}/gene_trees"
DNA_DIR="${OUTPUT_DIR}/dna"
PROTEIN_DIR="${OUTPUT_DIR}/protein"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Check if a PHYLIP alignment file contains duplicate sequences
# Returns 0 if duplicates found, 1 if no duplicates
check_duplicates() {
    local phylip_file="$1"

    if [ ! -f "$phylip_file" ]; then
        return 1  # File doesn't exist, no duplicates
    fi

    # Extract sequences from PHYLIP file (skip header line)
    # Sort them and check for duplicates
    local num_seqs=$(tail -n +2 "$phylip_file" | awk '{print $2}' | sort | uniq | wc -l)
    local total_seqs=$(tail -n +2 "$phylip_file" | wc -l)

    if [ "$num_seqs" -lt "$total_seqs" ]; then
        return 0  # Duplicates found
    else
        return 1  # No duplicates
    fi
}

# Count the number of leaves in a Newick tree
# Simple algorithm: count commas and add 1
count_tree_leaves() {
    local tree_file="$1"
    if [ ! -f "$tree_file" ]; then
        echo "0"
        return
    fi
    local tree_string=$(cat "$tree_file")
    local comma_count=$(echo "$tree_string" | tr -cd ',' | wc -c)
    echo $((comma_count + 1))
}

# ============================================================================
# FIND SOFTWARE BINARIES
# ============================================================================

# Find SimPhy
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -x "${SCRIPT_DIR}/SimPhy_1.0.2/bin/simphy_lnx64" ]; then
    SIMPHY="${SCRIPT_DIR}/SimPhy_1.0.2/bin/simphy_lnx64"
    echo "Found SimPhy at: ${SIMPHY}"
elif command -v simphy_lnx64 &> /dev/null; then
    SIMPHY="simphy_lnx64"
    echo "Found SimPhy in PATH"
else
    echo "Error: SimPhy (simphy_lnx64) not found!"
    echo "Please install SimPhy or add it to PATH"
    echo "Download from: https://github.com/adamallo/SimPhy"
    exit 1
fi

# Check for IQ-TREE 2 (AliSim)
if command -v iqtree2 &> /dev/null; then
    echo "Found IQ-TREE 2 (AliSim) in PATH"
else
    echo "Warning: IQ-TREE 2 (iqtree2) not found"
    echo "Phases 3 and 4 (sequence simulation) will fail"
    echo "Install from: http://www.iqtree.org/"
    echo "Or: sudo apt-get install iqtree (may need version 2.0+)"
fi

# Check for PhyML (only if ML estimation is enabled)
if [ ${RUN_ML_ESTIMATION} -eq 1 ]; then
    if ! command -v phyml &> /dev/null; then
        echo "Warning: PhyML not found in PATH"
        echo "Phases 3 and 4 (ML tree estimation) will fail"
        echo "Install with: sudo apt-get install phyml"
    fi
fi

echo ""

# ============================================================================
# SETUP
# ============================================================================

echo "=========================================="
echo "Simulation Pipeline Starting"
echo "=========================================="
echo ""
echo "Configuration:"
if [ ${INFER_ONLY} -eq 1 ]; then
    echo "  Mode: INFER-ONLY (only run PhyML on existing alignments)"
else
    echo "  Replicates: ${NUM_REPLICATES}"
    echo "  Leaf counts: ${NUM_LEAVES[@]}"
    echo "  Random seed: ${RANDOM_SEED}"
    if [ ${MIN_GENE_TREE_LEAVES} -gt 0 ]; then
        echo "  Gene tree filter: require exactly ${MIN_GENE_TREE_LEAVES} leaves"
    fi
fi
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Sequence types:"
if [ ${RUN_DNA} -eq 1 ]; then
    echo "    - DNA: YES"
else
    echo "    - DNA: NO"
fi
if [ ${RUN_PROTEIN} -eq 1 ]; then
    echo "    - Protein: YES"
else
    echo "    - Protein: NO"
fi
if [ ${INFER_ONLY} -eq 0 ]; then
    if [ ${RUN_ML_ESTIMATION} -eq 1 ]; then
        echo "  ML tree estimation: YES"
    else
        echo "  ML tree estimation: NO (alignments only)"
    fi
    echo "  Alignment lengths:"
    if [ ${RUN_DNA} -eq 1 ]; then
        echo "    - DNA: ${ALIGNMENT_LENGTH_DNA} bp"
    fi
    if [ ${RUN_PROTEIN} -eq 1 ]; then
        echo "    - Protein: ${ALIGNMENT_LENGTH_PROTEIN} aa"
    fi
    if [ ${ENABLE_INDELS} -eq 1 ]; then
        echo "  Indel simulation: YES"
        echo "    - Rates (ins,del): ${INDEL_RATES}"
        echo "    - Length distribution: ${INDEL_SIZE_DIST}"
    else
        echo "  Indel simulation: NO (substitutions only)"
    fi
fi
echo ""

# Create output directories
mkdir -p "${SPECIES_TREES_DIR}"
mkdir -p "${GENE_TREES_DIR}"
if [ ${RUN_DNA} -eq 1 ]; then
    mkdir -p "${DNA_DIR}"
fi
if [ ${RUN_PROTEIN} -eq 1 ]; then
    mkdir -p "${PROTEIN_DIR}"
fi

# ============================================================================
# PHASE 1: SIMULATE SPECIES TREES
# ============================================================================

echo "Phase 1: Simulating species trees..."
echo "--------------------------------------"

if [ ${INFER_ONLY} -eq 1 ]; then
    echo "  Skipped (infer-only mode)"
    echo ""
else

for num_leaves in "${NUM_LEAVES[@]}"; do
     echo "  Simulating ${NUM_REPLICATES} species trees with ${num_leaves} leaves..."

     echo "GENERATOR of species tree"

    "${SIMPHY}" -rs ${NUM_REPLICATES} -sb f:${SPECIATION_RATE} -sd f:${EXTINCTION_RATE} -sl f:${num_leaves} -stl f:${TREE_HEIGHT}  -sp f:100000 -o ${SPECIES_TREES_DIR}/lf${num_leaves}

    for i in $(seq -w 1 ${NUM_REPLICATES}); do
        echo REPL $i
        st=$(cat ${SPECIES_TREES_DIR}/lf${num_leaves}/$i/s_tree.trees)

        cat << nexus > ${SPECIES_TREES_DIR}/lf${num_leaves}_replicate$i.nex
#NEXUS
begin trees;
  tree tree1 = $st
end;
nexus

    done


    echo "----------DONE"

#     # Create R script for species tree simulation
#     cat > "${SPECIES_TREES_DIR}/simulate_species_trees_${num_leaves}.R" <<'RSCRIPT'
# library(TreeSim)

# # Parameters from command line or defaults
# args <- commandArgs(trailingOnly = TRUE)
# num_replicates <- as.integer(args[1])
# num_leaves <- as.integer(args[2])
# tree_height <- as.numeric(args[3])
# lambda <- as.numeric(args[4])  # speciation rate
# mu <- as.numeric(args[5])      # extinction rate
# output_dir <- args[6]
# random_seed <- as.integer(args[7])

# set.seed(random_seed)

# # Simulate species trees using general sampling approach
# for (i in 1:num_replicates) {
#     # Simulate birth-death tree
#     tree <- sim.bd.taxa.age(
#         n = num_leaves,
#         numbsim = 1,
#         lambda = lambda,
#         mu = mu,
#         age = tree_height,
#         mrca = TRUE
#     )[[1]]

#     # Save tree in Newick format
#     output_file <- file.path(output_dir,
#                             sprintf("species_tree_%d_leaves_%03d.nwk", num_leaves, i))
#     write.tree(tree, file = output_file)

#     # Save tree in Nexus format (required by SimPhy)
#     # SimPhy needs a simple Nexus format without TRANSLATE block
#     nexus_file <- file.path(output_dir,
#                            sprintf("species_tree_%d_leaves_%03d.nex", num_leaves, i))
    
#     tree$edge.length <- tree$edge.length 

#     # Remove root edge if it exists (SimPhy doesn't accept root branch lengths)
#     tree$root.edge <- NULL

#     newick_str <- write.tree(tree)

#     # SimPhy doesn't accept root branch lengths
#     # TreeSim may generate trees with root branches despite tree$root.edge <- NULL
#     # Remove root branch if present: (CONTENT):ROOT_BRANCH); -> (CONTENT);
#     # Pattern matches :NUMBER); at the very end and replaces with );
#     # newick_str <- sub(":[0-9.eE+-]+\\);$", ");", newick_str)

#     # Write simple Nexus format for SimPhy
#     cat("#NEXUS\n", file = nexus_file)
#     cat("begin trees;\n", file = nexus_file, append = TRUE)
#     cat(sprintf("  tree tree1 = %s\n", newick_str), file = nexus_file, append = TRUE)
#     cat("end;\n", file = nexus_file, append = TRUE)
# }

# cat(sprintf("Successfully simulated %d species trees with %d leaves\n",
#             num_replicates, num_leaves))
# RSCRIPT

#     # Run R script
#     Rscript "${SPECIES_TREES_DIR}/simulate_species_trees_${num_leaves}.R" \
#         "${NUM_REPLICATES}" \
#         "${num_leaves}" \
#         "${TREE_HEIGHT}" \
#         "${SPECIATION_RATE}" \
#         "${EXTINCTION_RATE}" \
#         "${SPECIES_TREES_DIR}" \
#         "${RANDOM_SEED}"

done




echo "Phase 1 completed!"
echo ""

fi  # End of INFER_ONLY check for Phase 1

# ============================================================================
# PHASE 2: SIMULATE GENE TREES
# ============================================================================

echo "Phase 2: Simulating gene trees with SimPhy..."
echo "--------------------------------------"

replicate_counter=0

for num_leaves in "${NUM_LEAVES[@]}"; do
    for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
        for pop_size in "${POPULATION_SIZES[@]}"; do

            config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
            config_dir="${GENE_TREES_DIR}/${config_name}"
            mkdir -p "${config_dir}"

            echo "  Configuration: ${config_name}"


            for i in $(seq -w 1 ${NUM_REPLICATES}); do
                species_tree="${SPECIES_TREES_DIR}/lf${num_leaves}_replicate$i.nex"

                output_folder="${config_dir}/replicate_$i"

                # Initialize gene tree filtering variables
                gene_tree_accepted=0
                gene_tree_retry=0

                # Loop until we get a gene tree with the required number of leaves (or reach max retries)
                while [ $gene_tree_accepted -eq 0 ]; do
                    # Generate unique seed for THIS specific simulation
                    # CRITICAL: Must vary by configuration AND replicate to get independent coalescent simulations
                    simphy_seed=$((RANDOM_SEED + replicate_counter + gene_tree_retry * 1000000))

                    # **Parameters:**
                    # - `-sr`: Species tree in Nexus format
                    # - `-rg 1`: Generate 1 gene tree per locus
                    # - `-rl F:1`: Generate 1 locus
                    # - `-si F:1`: Generate 1 individual per species
                    # - `-sp F:<pop_size>`: Effective population size (10⁷ or 5×10⁷)
                    # - `-su F:0.0004`: Substitution rate (scaled)
                    # - `-lb F:<scaled_dl_rate>`: Birth (duplication) rate (scaled)
                    # - `-ld F:lb`: Death (loss) rate = birth rate
                    # - `-ll 3`: Minimum number of lineages
                    # - `-hg LN:1.5,1`: Gene tree height follows lognormal distribution
                    # - `-oc 1`: Output coalescent trees
                    # - `-cs <unique_seed>`: Random seed (base_seed + replicate_number for independence)


                    while "${SIMPHY}" \
                        -sr "${species_tree}" \
                        -rg 1 \
                        -rl F:1 \
                        -si F:1 \
                        -sp "F:${pop_size}" \
                        -su F:0.0000000004 \
                        -lb "F:${dl_rate}" \
                        -ld F:lb \
                        -ll 3 \
                        -hg LN:1.5,1 \
                        -oc 1 \
                        -o "${output_folder}" \
                        -v 0 \
                        -cs ${simphy_seed}

                    do
                        break # hack, simphy returns undeterministic error sometimes
                    done

                    # Check if gene tree filtering is enabled
                    if [ ${MIN_GENE_TREE_LEAVES} -gt 0 ]; then
                        gene_tree_file="${output_folder}/1/g_trees1.trees"

                        if [ -f "$gene_tree_file" ]; then
                            gene_tree_leaves=$(count_tree_leaves "$gene_tree_file")

                            if [ "$gene_tree_leaves" -eq "${MIN_GENE_TREE_LEAVES}" ]; then
                                # Tree has correct number of leaves, accept it
                                gene_tree_accepted=1
                            else
                                # Tree has wrong number of leaves, retry
                                gene_tree_retry=$((gene_tree_retry + 1))

                                if [ $gene_tree_retry -ge ${MAX_GENE_TREE_RETRIES} ]; then
                                    echo "    Warning: replicate $i - max retries (${MAX_GENE_TREE_RETRIES}) reached. Gene tree has $gene_tree_leaves leaves (wanted ${MIN_GENE_TREE_LEAVES})"
                                    gene_tree_accepted=1  # Accept anyway to avoid infinite loop
                                fi
                            fi
                        else
                            echo "    Warning: Gene tree file not found: $gene_tree_file"
                            gene_tree_accepted=1  # Exit loop
                        fi
                    else
                        # No filtering, accept the tree
                        gene_tree_accepted=1
                    fi
                done

                if [ ${MIN_GENE_TREE_LEAVES} -gt 0 ] && [ $gene_tree_retry -gt 0 ]; then
                    gene_tree_file="${output_folder}/1/g_trees1.trees"
                    final_gene_tree_leaves=$(count_tree_leaves "$gene_tree_file")
                    if [ "$final_gene_tree_leaves" -eq "${MIN_GENE_TREE_LEAVES}" ]; then
                        echo "    Replicate $i: accepted after $gene_tree_retry retries ($final_gene_tree_leaves leaves)"
                    fi
                fi

                replicate_counter=$((replicate_counter + 1))
            done
        done
    done
done

echo "Phase 2 completed! Generated ${replicate_counter} gene trees."
echo ""



# ============================================================================
# PHASE 3: SIMULATE DNA SEQUENCES AND ESTIMATE GENE TREES
# ============================================================================

if [ ${RUN_DNA} -eq 1 ]; then

echo "Phase 3: Simulating DNA sequences and estimating gene trees..."
echo "--------------------------------------"

for num_leaves in "${NUM_LEAVES[@]}"; do
    for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
        for pop_size in "${POPULATION_SIZES[@]}"; do

            config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
            config_dir="${GENE_TREES_DIR}/${config_name}"
            dna_output_dir="${DNA_DIR}/${config_name}"
            mkdir -p "${dna_output_dir}"

            echo "  Processing DNA for: ${config_name}"

            for i in $(seq -w 1 ${NUM_REPLICATES}); do
                gene_tree_folder="${config_dir}/replicate_$i"
                dna_folder="${dna_output_dir}/replicate_$i"
                mkdir -p "${dna_folder}"

                # Read gene tree
                gene_tree_file="${gene_tree_folder}/1/g_trees1.trees"
                if [ -f "${gene_tree_file}" ]; then
                    gene_tree=$(cat "${gene_tree_file}")

                    # SimPhy with -su parameter already outputs branch lengths in substitutions per site
                    # No scaling needed - use gene tree directly
                    echo "${gene_tree}" > "${dna_folder}/gene_tree.nwk"

                    # Sample GTR parameters from Dirichlet distribution
                    # Returns: rate1 rate2 rate3 rate4 rate5 (rate6 normalized to 1.0)
                    gtr_rates=$(Rscript -e "
                        rates <- rgamma(6, c(12.776722, 20.869581, 5.647810, 9.863668, 30.679899, 3.199725), 1)
                        normalized <- rates / rates[6]
                        cat(paste(round(normalized[1:5], 6), collapse='/'))
                    ")

                    # Sample base frequencies from Dirichlet and normalize
                    base_freqs=$(Rscript -e "
                        f <- rgamma(4, c(113.48869, 69.02545, 78.66144, 99.83793), 1)
                        f <- f / sum(f)
                        cat(paste(round(f, 6), collapse='/'))
                    ")

                    # Sample gamma shape parameter from Lognormal
                    alpha=$(Rscript -e "cat(round(exp(rnorm(1, -0.470703916, 0.348667224)), 6))")

                    # Build AliSim model specification: GTR{rates}+F{freqs}+G4{alpha}
                    model="GTR{${gtr_rates}}+F{${base_freqs}}+G4{${alpha}}"

                    # Run AliSim to generate DNA sequences (with retry if duplicates found)
                    if command -v iqtree2 &> /dev/null; then
                        retry_count=0
                        success=0

                        while [ $retry_count -le $MAX_RETRIES ]; do
                            # Use different seed for each retry
                            # Force base-10 interpretation to avoid octal issues with leading zeros
                            current_seed=$((RANDOM_SEED + 10#$i + retry_count * 1000000))

                            # Build iqtree2 command with optional indel parameters
                            iqtree_cmd="iqtree2 --alisim ${dna_folder}/alignment -m ${model} -t ${dna_folder}/gene_tree.nwk --length ${ALIGNMENT_LENGTH_DNA}"
                            if [ ${ENABLE_INDELS} -eq 1 ]; then
                                # AliSim requires same distribution for insertion and deletion
                                iqtree_cmd="${iqtree_cmd} --indel ${INDEL_RATES} --indel-size ${INDEL_SIZE_DIST},${INDEL_SIZE_DIST}"
                            fi
                            iqtree_cmd="${iqtree_cmd} -seed ${current_seed} -af phy -quiet"

                            eval "${iqtree_cmd}" 2>/dev/null || true

                            # Rename output if it exists
                            if [ -f "${dna_folder}/alignment.phy" ]; then
                                mv "${dna_folder}/alignment.phy" "${dna_folder}/alignment_TRUE.phy"

                                # Check for duplicates if retry is enabled
                                if [ $MAX_RETRIES -gt 0 ] && check_duplicates "${dna_folder}/alignment_TRUE.phy"; then
                                    retry_count=$((retry_count + 1))
                                    if [ $retry_count -le $MAX_RETRIES ]; then
                                        echo "    DNA replicate ${i}: Duplicates found, retrying (attempt ${retry_count}/${MAX_RETRIES})"
                                        rm -f "${dna_folder}/alignment_TRUE.phy"
                                    else
                                        echo "    DNA replicate ${i}: Max retries reached, keeping alignment with duplicates"
                                        success=1
                                        break
                                    fi
                                else
                                    # No duplicates or retry disabled
                                    success=1
                                    break
                                fi
                            else
                                break  # AliSim failed
                            fi
                        done
                    else
                        echo "Warning: iqtree2 not found. Skipping sequence simulation for ${config_name} replicate ${i}"
                    fi
                fi

                # Estimate gene tree with PhyML
                if [ ${RUN_ML_ESTIMATION} -eq 1 ] && [ -f "${dna_folder}/alignment_TRUE.phy" ]; then
                    phyml -i "${dna_folder}/alignment_TRUE.phy" \
                          -m GTR \
                          -c 4 \
                          -a e \
                          -b 0 \
                          -o tlr \
                          --quiet 2>/dev/null

                    # Root the tree using midpoint rooting
                    if [ -f "${dna_folder}/alignment_TRUE.phy_phyml_tree.txt" ]; then
                        # This would use URec or a similar tool for midpoint-plateau rooting
                        # For now, we'll just copy the ML tree
                        cp "${dna_folder}/alignment_TRUE.phy_phyml_tree.txt" \
                           "${dna_folder}/ml_gene_tree.nwk"
                    fi
                fi
            done
        done
    done
done

echo "Phase 3 completed!"
echo ""

fi  # End of RUN_DNA check

# ============================================================================
# PHASE 4: SIMULATE PROTEIN SEQUENCES AND ESTIMATE GENE TREES
# ============================================================================

if [ ${RUN_PROTEIN} -eq 1 ]; then

echo "Phase 4: Simulating protein sequences and estimating gene trees..."
echo "--------------------------------------"

for num_leaves in "${NUM_LEAVES[@]}"; do
    for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
        for pop_size in "${POPULATION_SIZES[@]}"; do

            config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
            config_dir="${GENE_TREES_DIR}/${config_name}"
            protein_output_dir="${PROTEIN_DIR}/${config_name}"
            mkdir -p "${protein_output_dir}"

            echo "  Processing proteins for: ${config_name}"

            for i in $(seq -w 1 ${NUM_REPLICATES}); do
                gene_tree_folder="${config_dir}/replicate_$i"
                protein_folder="${protein_output_dir}/replicate_$i"
                mkdir -p "${protein_folder}"

                # Read gene tree
                gene_tree_file="${gene_tree_folder}/1/g_trees1.trees"
                if [ -f "${gene_tree_file}" ]; then
                    gene_tree=$(cat "${gene_tree_file}")

                    # SimPhy with -su parameter already outputs branch lengths in substitutions per site
                    # No scaling needed - use gene tree directly
                    echo "${gene_tree}" > "${protein_folder}/gene_tree.nwk"

                    # Sample gamma shape parameter from Lognormal
                    alpha=$(Rscript -e "cat(round(exp(rnorm(1, -0.470703916, 0.348667224)), 6))")

                    # Build AliSim model specification: WAG+G4{alpha}
                    model="WAG+G4{${alpha}}"

                    # Run AliSim to generate protein sequences (with retry if duplicates found)
                    if command -v iqtree2 &> /dev/null; then
                        retry_count=0
                        success=0

                        while [ $retry_count -le $MAX_RETRIES ]; do
                            # Use different seed for each retry
                            # Force base-10 interpretation to avoid octal issues with leading zeros
                            current_seed=$((RANDOM_SEED + 10#$i + retry_count * 1000000))

                            # Build iqtree2 command with optional indel parameters
                            iqtree_cmd="iqtree2 --alisim ${protein_folder}/alignment -m ${model} -t ${protein_folder}/gene_tree.nwk --length ${ALIGNMENT_LENGTH_PROTEIN}"
                            if [ ${ENABLE_INDELS} -eq 1 ]; then
                                # AliSim requires same distribution for insertion and deletion
                                iqtree_cmd="${iqtree_cmd} --indel ${INDEL_RATES} --indel-size ${INDEL_SIZE_DIST},${INDEL_SIZE_DIST}"
                            fi
                            iqtree_cmd="${iqtree_cmd} -seed ${current_seed} -af phy -quiet"

                            eval "${iqtree_cmd}" 2>/dev/null || true

                            # Rename output if it exists
                            if [ -f "${protein_folder}/alignment.phy" ]; then
                                mv "${protein_folder}/alignment.phy" "${protein_folder}/alignment_TRUE.phy"

                                # Check for duplicates if retry is enabled
                                if [ $MAX_RETRIES -gt 0 ] && check_duplicates "${protein_folder}/alignment_TRUE.phy"; then
                                    retry_count=$((retry_count + 1))
                                    if [ $retry_count -le $MAX_RETRIES ]; then
                                        echo "    Protein replicate ${i}: Duplicates found, retrying (attempt ${retry_count}/${MAX_RETRIES})"
                                        rm -f "${protein_folder}/alignment_TRUE.phy"
                                    else
                                        echo "    Protein replicate ${i}: Max retries reached, keeping alignment with duplicates"
                                        success=1
                                        break
                                    fi
                                else
                                    # No duplicates or retry disabled
                                    success=1
                                    break
                                fi
                            else
                                break  # AliSim failed
                            fi
                        done
                    else
                        echo "Warning: iqtree2 not found. Skipping sequence simulation for ${config_name} replicate ${i}"
                    fi
                fi

                # Estimate gene tree with PhyML (protein mode)
                if [ ${RUN_ML_ESTIMATION} -eq 1 ] && [ -f "${protein_folder}/alignment_TRUE.phy" ]; then
                    phyml -i "${protein_folder}/alignment_TRUE.phy" \
                          -d aa \
                          -m WAG \
                          -c 4 \
                          -a e \
                          -b 0 \
                          -o tlr \
                          --quiet 2>/dev/null

                    # Root the tree using midpoint rooting
                    if [ -f "${protein_folder}/alignment_TRUE.phy_phyml_tree.txt" ]; then
                        cp "${protein_folder}/alignment_TRUE.phy_phyml_tree.txt" \
                           "${protein_folder}/ml_gene_tree.nwk"
                    fi
                fi
            done
        done
    done
done

echo "Phase 4 completed!"
echo ""

fi  # End of RUN_PROTEIN check

# ============================================================================
# SUMMARY AND METADATA FILE
# ============================================================================

# Create summary file with all parameters and file statistics
SUMMARY_FILE="${OUTPUT_DIR}/simulation_summary.txt"

echo "Writing simulation summary to: ${SUMMARY_FILE}"

cat > "${SUMMARY_FILE}" << 'EOF_HEADER'
================================================================================
SIMULATION PIPELINE SUMMARY
================================================================================
EOF_HEADER

# Add timestamp
cat >> "${SUMMARY_FILE}" << EOF

Generated: $(date '+%Y-%m-%d %H:%M:%S')
Hostname: $(hostname)
Working Directory: $(pwd)

================================================================================
CONFIGURATION PARAMETERS
================================================================================

Pipeline Mode:
  Infer-only mode: $([ ${INFER_ONLY} -eq 1 ] && echo "YES (only ML tree inference)" || echo "NO (full simulation)")
  ML tree estimation: $([ ${RUN_ML_ESTIMATION} -eq 1 ] && echo "YES" || echo "NO")

Species Tree Parameters:
  Tree height: ${TREE_HEIGHT} years
  Speciation rate: ${SPECIATION_RATE} events/year
  Extinction rate: ${EXTINCTION_RATE} events/year
  Leaf set sizes: ${NUM_LEAVES[@]}

Gene Tree Parameters (SimPhy):
  Duplication/Loss rates: ${DUPLICATION_LOSS_RATES[@]} events/year
  Population sizes: ${POPULATION_SIZES[@]}
  Number of replicates: ${NUM_REPLICATES}
  Random seed (base): ${RANDOM_SEED}

Gene Tree Filtering:
  Minimum leaves required: $([ ${MIN_GENE_TREE_LEAVES} -gt 0 ] && echo "${MIN_GENE_TREE_LEAVES} (filtering enabled)" || echo "0 (no filtering)")
  Max retry attempts: ${MAX_GENE_TREE_RETRIES}

Sequence Simulation:
  Sequence types: $([ ${RUN_DNA} -eq 1 ] && echo -n "DNA " || echo -n "")$([ ${RUN_PROTEIN} -eq 1 ] && echo -n "PROTEIN" || echo -n "")
  DNA alignment length: ${ALIGNMENT_LENGTH_DNA} bp
  Protein alignment length: ${ALIGNMENT_LENGTH_PROTEIN} aa

Indel Simulation:
  Enabled: $([ ${ENABLE_INDELS} -eq 1 ] && echo "YES" || echo "NO")
EOF

if [ ${ENABLE_INDELS} -eq 1 ]; then
    cat >> "${SUMMARY_FILE}" << EOF
  Insertion/Deletion rates: ${INDEL_RATES}
  Length distribution: ${INDEL_SIZE_DIST}
EOF
fi

cat >> "${SUMMARY_FILE}" << EOF

Duplicate Sequence Handling:
  Max retry attempts: ${MAX_RETRIES} $([ ${MAX_RETRIES} -gt 0 ] && echo "(retry enabled)" || echo "(no retry)")

Output Directory: ${OUTPUT_DIR}

================================================================================
CONFIGURATION COMBINATIONS
================================================================================

Total configurations: $((${#NUM_LEAVES[@]} * ${#DUPLICATION_LOSS_RATES[@]} * ${#POPULATION_SIZES[@]}))

Configuration matrix:
EOF

for num_leaves in "${NUM_LEAVES[@]}"; do
    for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
        for pop_size in "${POPULATION_SIZES[@]}"; do
            config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
            echo "  - ${config_name}" >> "${SUMMARY_FILE}"
        done
    done
done

# Count generated files
if [ ${INFER_ONLY} -eq 0 ]; then
    cat >> "${SUMMARY_FILE}" << EOF

================================================================================
GENERATED FILES SUMMARY
================================================================================

Species Trees:
EOF

    species_tree_count=$(find "${SPECIES_TREES_DIR}" -name "*.nex" 2>/dev/null | wc -l)
    echo "  Total species trees: ${species_tree_count}" >> "${SUMMARY_FILE}"

    cat >> "${SUMMARY_FILE}" << EOF

Gene Trees:
EOF

    gene_tree_count=$(find "${GENE_TREES_DIR}" -name "g_trees1.trees" 2>/dev/null | wc -l)
    echo "  Total gene trees: ${gene_tree_count}" >> "${SUMMARY_FILE}"
fi

if [ ${RUN_DNA} -eq 1 ]; then
    cat >> "${SUMMARY_FILE}" << EOF

DNA Sequences and Trees:
EOF

    dna_alignment_count=$(find "${DNA_DIR}" -name "alignment_TRUE.phy" 2>/dev/null | wc -l)
    dna_ml_tree_count=$(find "${DNA_DIR}" -name "ml_gene_tree.nwk" 2>/dev/null | wc -l)
    echo "  DNA alignments: ${dna_alignment_count}" >> "${SUMMARY_FILE}"
    if [ ${RUN_ML_ESTIMATION} -eq 1 ]; then
        echo "  DNA ML trees: ${dna_ml_tree_count}" >> "${SUMMARY_FILE}"
    fi
fi

if [ ${RUN_PROTEIN} -eq 1 ]; then
    cat >> "${SUMMARY_FILE}" << EOF

Protein Sequences and Trees:
EOF

    protein_alignment_count=$(find "${PROTEIN_DIR}" -name "alignment_TRUE.phy" 2>/dev/null | wc -l)
    protein_ml_tree_count=$(find "${PROTEIN_DIR}" -name "ml_gene_tree.nwk" 2>/dev/null | wc -l)
    echo "  Protein alignments: ${protein_alignment_count}" >> "${SUMMARY_FILE}"
    if [ ${RUN_ML_ESTIMATION} -eq 1 ]; then
        echo "  Protein ML trees: ${protein_ml_tree_count}" >> "${SUMMARY_FILE}"
    fi
fi

# Software versions
cat >> "${SUMMARY_FILE}" << EOF

================================================================================
SOFTWARE VERSIONS
================================================================================

EOF

if command -v simphy_lnx64 &> /dev/null; then
    echo "SimPhy: $(simphy_lnx64 -h 2>&1 | grep -i "version" | head -1 || echo "Version information not available")" >> "${SUMMARY_FILE}"
else
    echo "SimPhy: Not found in PATH" >> "${SUMMARY_FILE}"
fi

if command -v iqtree2 &> /dev/null; then
    iqtree_version=$(iqtree2 --version 2>&1 | grep -i "version" | head -1 || echo "Unknown")
    echo "IQ-TREE 2: ${iqtree_version}" >> "${SUMMARY_FILE}"
else
    echo "IQ-TREE 2: Not found in PATH" >> "${SUMMARY_FILE}"
fi

if command -v phyml &> /dev/null; then
    phyml_version=$(phyml --version 2>&1 | head -1 || echo "Unknown")
    echo "PhyML: ${phyml_version}" >> "${SUMMARY_FILE}"
else
    echo "PhyML: Not found in PATH" >> "${SUMMARY_FILE}"
fi

cat >> "${SUMMARY_FILE}" << EOF

================================================================================
DIRECTORY STRUCTURE
================================================================================

${OUTPUT_DIR}/
├── simulation_summary.txt (this file)
├── species_trees/
EOF

if [ ${INFER_ONLY} -eq 0 ]; then
    for num_leaves in "${NUM_LEAVES[@]}"; do
        echo "│   ├── lf${num_leaves}_replicate*.nex" >> "${SUMMARY_FILE}"
    done
fi

cat >> "${SUMMARY_FILE}" << EOF
├── gene_trees/
EOF

for num_leaves in "${NUM_LEAVES[@]}"; do
    for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
        for pop_size in "${POPULATION_SIZES[@]}"; do
            config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
            echo "│   ├── ${config_name}/" >> "${SUMMARY_FILE}"
            echo "│   │   └── replicate_*/1/g_trees1.trees" >> "${SUMMARY_FILE}"
        done
    done
done

if [ ${RUN_DNA} -eq 1 ]; then
    cat >> "${SUMMARY_FILE}" << EOF
├── dna/
EOF
    for num_leaves in "${NUM_LEAVES[@]}"; do
        for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
            for pop_size in "${POPULATION_SIZES[@]}"; do
                config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
                echo "│   ├── ${config_name}/" >> "${SUMMARY_FILE}"
                echo "│   │   └── replicate_*/{alignment_TRUE.phy, ml_gene_tree.nwk}" >> "${SUMMARY_FILE}"
            done
        done
    done
fi

if [ ${RUN_PROTEIN} -eq 1 ]; then
    cat >> "${SUMMARY_FILE}" << EOF
└── protein/
EOF
    for num_leaves in "${NUM_LEAVES[@]}"; do
        for dl_rate in "${DUPLICATION_LOSS_RATES[@]}"; do
            for pop_size in "${POPULATION_SIZES[@]}"; do
                config_name="leaves${num_leaves}_dl${dl_rate}_ps${pop_size}"
                echo "    ├── ${config_name}/" >> "${SUMMARY_FILE}"
                echo "    │   └── replicate_*/{alignment_TRUE.phy, ml_gene_tree.nwk}" >> "${SUMMARY_FILE}"
            done
        done
    done
fi

cat >> "${SUMMARY_FILE}" << EOF

================================================================================
NOTES
================================================================================

- All parameters are stored in this summary file for reproducibility
- Random seeds vary by replicate: base_seed + replicate_number
- Gene tree branch lengths are in substitutions per site (output by SimPhy)
- ML trees are unrooted (PhyML output)

For detailed methodology, see PIPELINE.md

================================================================================
EOF

echo "=========================================="
echo "Pipeline Completed Successfully!"
echo "=========================================="
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "Summary file: ${SUMMARY_FILE}"
echo ""
echo "Generated data:"
echo "  - Species trees: ${SPECIES_TREES_DIR}"
echo "  - Gene trees: ${GENE_TREES_DIR}"
echo "  - DNA sequences and ML trees: ${DNA_DIR}"
echo "  - Protein sequences and ML trees: ${PROTEIN_DIR}"
echo ""
echo "Configuration summary:"
echo "  - Replicates per configuration: ${NUM_REPLICATES}"
echo "  - Leaf set sizes: ${NUM_LEAVES[@]}"
echo "  - DL rates: ${DUPLICATION_LOSS_RATES[@]}"
echo "  - Population sizes: ${POPULATION_SIZES[@]}"
echo "  - Total configurations: $((${#NUM_LEAVES[@]} * ${#DUPLICATION_LOSS_RATES[@]} * ${#POPULATION_SIZES[@]}))"
echo ""
