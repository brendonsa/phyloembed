# Simulation Pipeline for MSAs and Gene Trees

This pipeline simulates multiple sequence alignments (MSAs) and gene trees for both DNA and protein sequences, following the methodology described in Molloy and Warnow (2020).

## Overview

The pipeline consists of four phases:

1. **Species Tree Simulation** - Generate ultrametric species trees using SimPhy's birth-death process
2. **Gene Tree Simulation** - Simulate gene trees with duplication/loss and incomplete lineage sorting (ILS) using SimPhy
3. **DNA Sequence Simulation** - Generate DNA sequences, MSAs, and estimate gene trees using maximum likelihood (AliSim + PhyML)
4. **Protein Sequence Simulation** - Generate protein sequences, MSAs, and estimate gene trees using maximum likelihood (AliSim + PhyML)

## Requirements

### Software Dependencies

- **SimPhy** (version 1.0.2+)
  - Download from: https://github.com/adamallo/SimPhy
  - Used for both species tree and gene tree simulation

- **IQ-TREE 2** (version 2.2.0+, for AliSim sequence simulation)
  - **IMPORTANT**: Repository versions may be too old (e.g., Ubuntu apt has 2.0.7). AliSim requires 2.2.0+
  - Download latest from: http://www.iqtree.org/ or https://github.com/iqtree/iqtree2/releases
  - Manual install recommended:
    ```bash
    wget https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz
    tar -xzf iqtree-2.3.6-Linux-intel.tar.gz
    sudo cp iqtree-2.3.6-Linux-intel/bin/iqtree2 /usr/local/bin/
    ```

- **PhyML** (version 3.1+)
  - Download from: http://www.atgc-montpellier.fr/phyml/
  - Ubuntu/Debian: `sudo apt-get install phyml`

- **URec** (optional, for midpoint-plateau rooting)
  - Download from: http://www.lcqb.upmc.fr/URec/

## Pipeline Parameters

### Species Tree Parameters
- **Tree height**: 1,800,000,337.5 years
- **Speciation rate**: 1.8 × 10⁻⁹ events/year
- **Extinction rate**: 0 events/year
- **Number of leaves**: 12 or 20

### Gene Tree Parameters (SimPhy)
- **Duplication/Loss rates**: {10⁻¹⁰, 2×10⁻¹⁰, 5×10⁻¹⁰}
- **Population sizes**: {10⁷, 5×10⁷} (low and medium ILS)
- **Substitution rate**: 4 × 10⁻¹⁰ substitutions/site/year

### Sequence Simulation Parameters

#### DNA Sequences (AliSim/IQ-TREE 2)
- **Model**: GTR+Γ
- **Substitution rates (AC, AG, AT, CG, CT, GT)**: Sampled from Dirichlet(12.776722, 20.869581, 5.647810, 9.863668, 30.679899, 3.199725)
- **Nucleotide frequencies (T, C, A, G)**: Sampled from Dirichlet(113.48869, 69.02545, 78.66144, 99.83793)
- **Alpha (Γ shape)**: Sampled from Lognormal(-0.470703916, 0.348667224)
- **Alignment length**: 1,000 bp (default, configurable via `--dna-length`)

#### Protein Sequences (AliSim/IQ-TREE 2)
- **Model**: WAG+Γ
- **Alpha (Γ shape)**: Sampled from Lognormal(-0.470703916, 0.348667224)
- **Alignment length**: 333 amino acids (default, configurable via `--protein-length`)

### Gene Tree Estimation
- **DNA**: PhyML with GTR+Γ model
- **Protein**: PhyML with WAG+Γ model
- **Rooting**: Midpoint-plateau rooting (URec)

## Usage

### Running the Full Pipeline

#### Bash Script

```bash
# Run with default parameters (100 replicates, 12 and 20 leaves, both DNA and protein)
./simulate_pipeline.sh

# Quick test: 10 replicates, 12 leaves only, both types
./simulate_pipeline.sh -n 10 -l 12

# Run 50 replicates with DNA only
./simulate_pipeline.sh -n 50 -t dna

# Run with 12, 20, and 50 leaves, protein only
./simulate_pipeline.sh -l 12,20,50 -t protein

# Generate alignments only (skip ML tree estimation, ~3x faster)
./simulate_pipeline.sh -n 100 -l 12 -m

# Later: infer ML trees from existing alignments
./simulate_pipeline.sh -o simulation_output -n 100 -l 12 -i

# Custom alignment lengths
./simulate_pipeline.sh -n 100 -l 12 --protein-length 500
./simulate_pipeline.sh -n 100 -l 12 --dna-length 2000 --protein-length 500

# With indel simulation
./simulate_pipeline.sh -n 100 -l 12 --indel-model realistic
./simulate_pipeline.sh -n 100 -l 12 --protein-length 500 --indel-model realistic

# Show help
./simulate_pipeline.sh -h
```

**Available Options:**
- `-n NUM` - Number of replicates per configuration (default: 100)
- `-l LEAVES` - Comma-separated list of leaf counts (default: 12,20)
- `-t TYPE` - Sequence type: dna, protein, or both (default: both)
- `-o DIR` - Output directory (default: simulation_output)
- `-s SEED` - Random seed (default: 22)
- `-r MAX` - Retry if duplicate sequences found, up to MAX attempts (default: 0)
- `-f NUM` - Filter gene trees: require exactly NUM leaves (default: 0 = no filtering)
- `-m` - Skip ML tree estimation (generate alignments only, ~3x faster)
- `-i` - Infer-only mode: run only PhyML on existing alignments
- `--dna-length NUM` - DNA alignment length in base pairs (default: 1000)
- `--protein-length NUM` - Protein alignment length in amino acids (default: 333)
- `--indel-model MODEL` - Indel preset: noindels, realistic, conservative, or highrate
- `--indel INS,DEL` - Enable indel simulation with specified rates (e.g., 0.03,0.09)
- `--indel-size DIST` - Indel length distribution (e.g., "POW{1.7/50}")
- `-h` - Show help message

#### Python Script

```bash
# Run with default parameters
./simulate_pipeline.py

# Quick test: 10 replicates, 12 leaves only
./simulate_pipeline.py --replicates 10 --leaves 12

# Run 50 replicates with DNA only
./simulate_pipeline.py -n 50 -t dna

# Run with multiple leaf counts, protein only
./simulate_pipeline.py -l 12 20 50 -t protein

# Show help
./simulate_pipeline.py -h
```

**Available Options:**
- `-n, --replicates NUM` - Number of replicates per configuration (default: 100)
- `-l, --leaves LEAVES` - Leaf counts (space-separated, default: 12 20)
- `-t, --type TYPE` - Sequence type: dna, protein, or both (default: both)
- `-o, --output DIR` - Output directory (default: simulation_output)
- `-s, --seed SEED` - Random seed (default: 22)
- `-h, --help` - Show help message

### Two-Stage Workflow (Alignments First, ML Trees Later)

For faster testing or when you want to defer ML tree inference, you can split the pipeline into two stages:

**Stage 1: Generate alignments only** (using `-m` flag)
```bash
# Fast: generate only alignments, skip PhyML (~3x faster)
./simulate_pipeline.sh -n 100 -l 12,20 -o my_output -m
```

This creates:
- Species trees
- Gene trees (with duplication/loss and ILS)
- DNA and protein alignments
- **No ML trees** (PhyML skipped)

**Stage 2: Infer ML trees from existing alignments** (using `-i` flag)
```bash
# Later: run only PhyML on existing alignments
./simulate_pipeline.sh -o my_output -n 100 -l 12,20 -i
```

This:
- **Skips** all simulation phases (species trees, gene trees, sequences)
- **Runs only PhyML** on existing `alignment_TRUE.phy` files
- Generates `ml_gene_tree.nwk` and PhyML statistics

**Important:** When using `-i`, you must specify the **same parameters** (`-n`, `-l`) used in Stage 1.

**Use Cases:**
- Quick testing: Generate alignments fast, infer trees later if needed
- Batch processing: Generate many alignments, then distribute ML inference across machines
- Parameter exploration: Generate alignments once, experiment with different ML settings

### Gene Tree Filtering by Leaf Count

Due to duplication and loss events in Phase 2, gene trees may have varying numbers of leaves. The `-f NUM` option allows you to filter gene trees to have exactly NUM leaves.

**How it works:**
```bash
# Require all gene trees to have exactly 10 leaves
./simulate_pipeline.sh -n 100 -l 12 -f 10
```

When `-f NUM` is specified:
- SimPhy generates a gene tree
- The pipeline counts the number of leaves
- If the tree has exactly NUM leaves, it's accepted
- If not, SimPhy regenerates with a different seed (up to 1000 attempts)
- Progress messages show how many retries were needed

**Use Cases:**
- **Consistent tree sizes**: Ensure all replicates have the same number of taxa for fair comparison
- **Filtering DL effects**: Remove extreme duplication/loss outcomes
- **Controlled experiments**: Study specific tree sizes while maintaining DL and ILS processes

**Example output:**
```
Configuration: leaves12_dl1e-10_ps1e7
Replicate 1: accepted after 15 retries (10 leaves)
Replicate 2: accepted after 3 retries (10 leaves)
```

**Note:** If a tree with exactly NUM leaves cannot be generated after 1000 attempts, the last tree is accepted with a warning message.

### Customizing Alignment Lengths

Alignment lengths can be customized via command-line options:

```bash
# Custom protein alignment length (e.g., for longer protein domains)
./simulate_pipeline.sh -n 100 -l 12 --protein-length 500

# Custom DNA alignment length
./simulate_pipeline.sh -n 100 -l 12 --dna-length 2000

# Both custom lengths
./simulate_pipeline.sh -n 100 -l 12 --dna-length 1500 --protein-length 500
```

**Effects on simulation:**
- **Without indels**: Final alignment length equals the specified value exactly
- **With indels**: Final alignment length varies due to insertions and deletions
  - Example: 1000 bp with `--indel 0.03,0.09` → ~940 bp on average

### Customizing Other Advanced Parameters

For duplication/loss rates, population sizes, and other parameters, edit the configuration section at the top of the script:

```bash
# In simulate_pipeline.sh:
DUPLICATION_LOSS_RATES=(1e-10 2e-10 5e-10)
POPULATION_SIZES=(1e7 5e7)
```

Or use the `config.example.sh` file:

```bash
cp config.example.sh config.sh
# Edit config.sh with your parameters
# Then source it in your script or use as reference
```

## Output Structure

```
simulation_output/
├── simulation_summary.txt              # Automatic summary of all parameters and files
├── species_trees/
│   ├── species_tree_12_leaves_001.nwk
│   ├── species_tree_12_leaves_001.nex
│   └── ...
├── gene_trees/
│   ├── leaves12_dl1e-10_ps1e7/
│   │   ├── replicate_001/
│   │   └── ...
│   └── ...
├── dna/
│   ├── leaves12_dl1e-10_ps1e7/
│   │   ├── replicate_001/
│   │   │   ├── alignment_TRUE.phy
│   │   │   ├── alignment.unaligned.fa
│   │   │   ├── alignment_TRUE.phy_phyml_tree.txt
│   │   │   └── ml_gene_tree.nwk
│   │   └── ...
│   └── ...
└── protein/
    ├── leaves12_dl1e-10_ps1e7/
    │   ├── replicate_001/
    │   │   ├── alignment_TRUE.phy
    │   │   ├── alignment.unaligned.fa
    │   │   ├── alignment_TRUE.phy_phyml_tree.txt
    │   │   └── ml_gene_tree.nwk
    │   └── ...
    └── ...
```

### Simulation Summary File

Each pipeline run automatically generates `simulation_summary.txt` in the output directory. This file contains:

- **Complete parameter documentation**: All command-line options, species tree parameters, gene tree settings, sequence lengths, indel parameters, etc.
- **File statistics**: Counts of all generated species trees, gene trees, alignments, and ML trees
- **Software versions**: SimPhy, IQ-TREE 2, and PhyML versions used
- **Configuration matrix**: List of all parameter combinations
- **Timestamp and environment**: When and where the simulation was run

This summary file is essential for:
- **Reproducibility**: Contains all information needed to reproduce the simulation
- **Documentation**: Quick reference for what parameters were used
- **Sharing**: Easy to share simulation details with collaborators
- **Verification**: Confirm all expected files were generated

**Example excerpt:**
```
================================================================================
CONFIGURATION PARAMETERS
================================================================================

Species Tree Parameters:
  Tree height: 1800000337.5 years
  Speciation rate: 1.8e-9 events/year
  Leaf set sizes: 12

Gene Tree Parameters (SimPhy):
  Duplication/Loss rates: 1e-10 2e-10 5e-10 events/year
  Population sizes: 1e7 5e7
  Number of replicates: 100

Sequence Simulation:
  DNA alignment length: 1000 bp
  Protein alignment length: 333 aa

Indel Simulation:
  Enabled: YES
  Insertion/Deletion rates: 0.03,0.09
  Length distribution: POW{1.7/50}

================================================================================
GENERATED FILES SUMMARY
================================================================================

Species Trees:
  Total species trees: 1

Gene Trees:
  Total gene trees: 6

Protein Sequences and Trees:
  Protein alignments: 6
  Protein ML trees: 6
```

## Configuration Combinations

With default parameters, the pipeline generates:
- 2 leaf set sizes (12, 20)
- 3 duplication/loss rates
- 2 population sizes
- **Total: 12 different configurations**
- Each with 100 replicates
- **Total datasets: 1,200 per sequence type (DNA/protein)**

## File Descriptions

### Summary File
- `simulation_summary.txt` - Automatically generated file containing all parameters, file counts, software versions, and configuration details for reproducibility

### Species Trees
- `.nwk` - Newick format (for general use)
- `.nex` - Nexus format (required by SimPhy)

### Gene Trees
- SimPhy output in `gene_trees/` directory
- Contains true gene trees with duplication, loss, and ILS events

### DNA/Protein Sequences
- `alignment_TRUE.phy` - True alignment (Phylip format)
- `alignment.unaligned.fa` - Unaligned sequences in FASTA format (automatically generated by AliSim, all gaps removed)
- `alignment_TRUE.phy_phyml_tree.txt` - ML tree estimated by PhyML
- `ml_gene_tree.nwk` - Rooted ML gene tree

## Biological Parameters

The parameters are based on fungal dataset from Rasmussen and Kellis (2011), specifically *Saccharomycetales* fungi, with some elevated parameters to test algorithmic performance in challenging conditions.

### Parameter Levels

**Duplication/Loss:**
- Low: 10⁻¹⁰ (biological data)
- Medium: 2×10⁻¹⁰
- High: 5×10⁻¹⁰

**Incomplete Lineage Sorting (ILS):**
- Low: Population size = 10⁷
- Medium: Population size = 5×10⁷

## Notes

1. **Runtime**: Full pipeline with default parameters may take several hours depending on hardware
2. **Disk Space**: Ensure adequate disk space (several GB for 1,200+ datasets)
3. **Random Seed**: Set to 22 by default for reproducibility. Each replicate uses seed+replicate_number for independent stochastic events
4. **Error Handling**: Script exits on first error (`set -e`)
5. **Critical Fixes (2025-12-23)**: Fixed duplication/loss rate scaling and random seed propagation. See `BUGFIX_DL_AND_SEED.md` for details. **Re-run old simulations if generated before this date.**

## Troubleshooting

### Common Issues

1. **SimPhy not found**
   - Ensure `simphy_lnx64` is in your PATH or provide full path in script
   - Local installation should be in `SimPhy_1.0.2/bin/`

2. **IQ-TREE 2 not found**
   - Install: `sudo apt-get install iqtree` (Ubuntu/Debian)
   - Or download from: http://www.iqtree.org/
   - Ensure version 2.0+ for AliSim support

3. **PhyML not found**
   - Install PhyML: `sudo apt-get install phyml` (Ubuntu/Debian)

**For detailed troubleshooting, see `KNOWN_ISSUES.md`**

## References

- Molloy, E.K. and Warnow, T. (2020). TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees.
- Rasmussen, M.D. and Kellis, M. (2011). A Bayesian approach for fast and accurate gene tree reconstruction.
- SimPhy: Mallo, D., De Oliveira Martins, L., and Posada, D. (2016). SimPhy: Phylogenomic simulation of gene, locus, and species trees. Systematic Biology, 65(2):334-344.
- **AliSim**: Ly-Trong, N. et al. (2022). AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era. Molecular Biology and Evolution, 39(5):msac092.
- **IQ-TREE 2**: Minh, B.Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5):1530-1534.
- PhyML: Guindon, S. et al. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology, 59(3):307-321.

## License

This pipeline is provided as-is for research purposes.
