# Simulation Pipeline Description

Our simulation procedure consists of four major phases: (i) simulating species trees, (ii) simulating gene trees with duplication/loss and incomplete lineage sorting, (iii) simulating DNA sequences and estimating gene trees, and (iv) simulating protein sequences and estimating gene trees. The selection of parameters in all phases is based on the simulation study by Molloy and Warnow (2020), which uses parameters derived from a fungal dataset presented by Rasmussen and Kellis (2011).

## Phase 1: Simulating Species Trees

Species trees are simulated using **SimPhy** version 1.0.2+ with a birth-death process using parameters from Molloy and Warnow (2020):

- **Tree height**: 1,800,000,337.5 years
- **Speciation rate**: λ = 1.8 × 10⁻⁹ events/year
- **Extinction rate**: μ = 0 events/year
- **Number of leaves**: 12 or 20 (configurable via command-line)

SimPhy produces **ultrametric species trees**, where all leaves are equidistant from the root.

### SimPhy Command for Species Trees

```bash
simphy_lnx64 \
  -rs <num_replicates> \
  -sb f:<speciation_rate> \
  -sd f:<extinction_rate> \
  -sl f:<num_leaves> \
  -stl f:<tree_height> \
  -sp f:100000 \
  -o <output_dir>
```

Trees are extracted from SimPhy output (`s_tree.trees`) and saved in simple Nexus (.nex) format for use in Phase 2.

## Phase 2: Simulating Gene Trees

Gene trees are simulated using **SimPhy** version 1.0.2+ to model:
- **Duplication and loss (DL) events** along species tree branches
- **Incomplete lineage sorting (ILS)** through coalescent simulation

### Key Parameters

- **Substitution rate**: 4 × 10⁻¹⁰ substitutions/site/year
- **Duplication/loss rates**: {10⁻¹⁰, 2×10⁻¹⁰, 5×10⁻¹⁰} events/year
- **Population sizes**: {10⁷, 5×10⁷}

SimPhy uses these rates directly (in years) and outputs gene tree branch lengths already scaled to **substitutions per site**.

### SimPhy Command

```bash
simphy_lnx64 \
  -sr <species_tree.nex> \
  -rg 1 \
  -rl F:1 \
  -si F:1 \
  -sp F:<population_size> \
  -su F:0.0000000004 \
  -lb F:<dl_rate> \
  -ld F:lb \
  -ll 3 \
  -hg LN:1.5,1 \
  -oc 1 \
  -o <output_folder> \
  -v 0 \
  -cs <unique_seed>
```

**Parameters:**
- `-sr`: Species tree in Nexus format
- `-rg 1`: Generate 1 gene tree per locus
- `-rl F:1`: Generate 1 locus
- `-si F:1`: Generate 1 individual per species
- `-sp F:<pop_size>`: Effective population size (10⁷ or 5×10⁷)
- `-su F:0.0000000004`: Substitution rate (4×10⁻¹⁰ substitutions/site/year)
- `-lb F:<dl_rate>`: Birth (duplication) rate (e.g., 1e-10, 2e-10, or 5e-10 events/year)
- `-ld F:lb`: Death (loss) rate = birth rate
- `-ll 3`: Minimum number of lineages
- `-hg LN:1.5,1`: Gene tree height follows lognormal distribution
- `-oc 1`: Output coalescent trees
- `-cs <unique_seed>`: Random seed (varies by configuration and replicate for independence)

### Parameter Combinations

Following Molloy and Warnow (2020), we use:

**Duplication/Loss Rates** (3 levels):
- Low: 10⁻¹⁰ events/year (biological data from *Saccharomycetales* fungi)
- Medium: 2×10⁻¹⁰ events/year
- High: 5×10⁻¹⁰ events/year

**Population Sizes** (2 levels):
- Low ILS: 10⁷ (corresponding to low incomplete lineage sorting)
- Medium ILS: 5×10⁷ (corresponding to medium incomplete lineage sorting)

This yields **6 parameter combinations** per leaf-set size, and **12 total configurations** for 2 leaf-set sizes.

### Critical Implementation Details

1. **Direct rate usage**: SimPhy uses rates directly in years (e.g., 1e-10, not scaled)
2. **Independent seeds**: Each configuration and replicate uses a unique seed to ensure independent stochastic events
3. **Branch lengths in substitutions**: SimPhy outputs gene tree branch lengths in **substitutions per site** (not time units)
4. **Output**: SimPhy produces three tree files per simulation:
   - `s_tree.trees`: Species tree
   - `l_trees.trees`: Locus tree (after duplication/loss)
   - `g_trees1.trees`: Gene tree (after ILS coalescence)

## Phase 3: Simulating DNA Sequences and Estimating Gene Trees

DNA sequences are simulated using **AliSim** (part of IQ-TREE 2) with parameters from Molloy and Warnow (2020).

Gene trees from SimPhy (Phase 2) are used directly - no branch length conversion is needed since SimPhy outputs branch lengths in **substitutions per site**.

### DNA Sequence Parameters

**Model**: GTR+Γ (General Time Reversible with Gamma rate heterogeneity)

**Optional: Insertion-Deletion (Indel) Model**

The pipeline supports realistic indel simulation using empirically-derived parameters (Ly-Trong et al. 2022). Indels can be enabled with the `--indel` flag.

**Substitution Rates** (AC, AG, AT, CG, CT, GT):
- Sampled from Dirichlet(12.776722, 20.869581, 5.647810, 9.863668, 30.679899, 3.199725)
- Independent sample for each gene tree

**Nucleotide Frequencies** (T, C, A, G):
- Sampled from Dirichlet(113.48869, 69.02545, 78.66144, 99.83793)
- Normalized to sum to 1.0
- Independent sample for each gene tree

**Gamma Shape Parameter** (α):
- Sampled from Lognormal(-0.470703916, 0.348667224)
- Controls rate heterogeneity across sites

**Alignment Length**: 1000 bp (default, configurable via `--dna-length`)

### AliSim Command Example

```bash
# Sample parameters
gtr_rates="2.345678/3.456789/0.876543/1.234567/4.567890"  # Relative to GT=1.0
base_freqs="0.301822/0.211993/0.225259/0.260926"         # Sum to 1.0
alpha="0.879134"

# Build model and run
model="GTR{${gtr_rates}}+F{${base_freqs}}+G4{${alpha}}"

# Without indels (default)
iqtree2 --alisim alignment \
  -m "${model}" \
  -t gene_tree.nwk \
  --length 1000 \
  -seed 22 \
  -af phy

# With indels (empirical parameters)
iqtree2 --alisim alignment \
  -m "${model}" \
  -t gene_tree.nwk \
  --length 1000 \
  --indel 0.03,0.09 \
  --indel-size POW 1.7 50 \
  -seed 22 \
  -af phy
```

### Maximum Likelihood Gene Tree Estimation

Gene trees are estimated from true alignments using **PhyML** v3.1+ with:

**Model**: GTR+Γ
```bash
phyml -i alignment_TRUE.phy \
  -m GTR \
  -c 4 \
  -a e \
  -b 0 \
  -o tlr \
  --quiet
```

**Parameters:**
- `-m GTR`: General Time Reversible model
- `-c 4`: 4 discrete gamma rate categories
- `-a e`: Estimate gamma shape parameter
- `-b 0`: No bootstrap
- `-o tlr`: Optimize topology, branch lengths, and substitution rates

**Output**: Unrooted ML tree in `alignment_TRUE.phy_phyml_tree.txt`

**Note**: ML tree estimation can be skipped using the `-m` flag (alignments only), or run separately using the `-i` flag (infer-only mode). See "Two-Stage Workflow" section below for details.

### Expected ML Tree Quality

For 1000 bp alignments with appropriate parameters:
- **Tree size**: 5-15 substitutions per site
- **Parsimony score**: 1000-2500 substitutions
- **Log-likelihood**: Typically -5000 to -7000
- **Branch lengths**: Reasonable range (0.001-5.0)

## Phase 4: Simulating Protein Sequences and Estimating Gene Trees

Protein sequences are simulated similarly to DNA, but using amino acid substitution models.

### Protein Sequence Parameters

**Model**: WAG+Γ (Whelan and Goldman model with Gamma rate heterogeneity)

**Gamma Shape Parameter** (α):
- Sampled from Lognormal(-0.470703916, 0.348667224)
- Same distribution as DNA, independent samples

**Alignment Length**: 333 amino acids (default, configurable via `--protein-length`)
- Default value chosen as ≈ 1000 bp / 3 to match DNA length
- Can be set independently from DNA length

### AliSim Command Example

```bash
# Sample gamma shape parameter
alpha="0.625"

# Build model and run
model="WAG+G4{${alpha}}"

iqtree2 --alisim alignment \
  -m "${model}" \
  -t gene_tree.nwk \
  --length 333 \
  -seed 22 \
  --out-format phylip
```

### Maximum Likelihood Estimation

Protein trees estimated using **PhyML** with:

```bash
phyml -i alignment_TRUE.phy \
  -d aa \
  -m WAG \
  -c 4 \
  -a e \
  -b 0 \
  -o tlr \
  --quiet
```

**Parameters:**
- `-d aa`: Amino acid sequences
- `-m WAG`: WAG substitution matrix
- Other parameters same as DNA

**Note**: ML tree estimation can be skipped using the `-m` flag (alignments only), or run separately using the `-i` flag (infer-only mode). See "Two-Stage Workflow" section below for details.

## Insertion-Deletion (Indel) Modeling

### Overview

The pipeline optionally supports insertion-deletion (indel) simulation using **AliSim** (Ly-Trong et al. 2022). While substitution-only models (GTR, WAG) are simpler and computationally faster, indel modeling provides more biologically realistic sequence evolution and is important for testing phylogenetic methods under realistic conditions.

### When to Use Indel Modeling

**Use indel simulation when:**
- Testing alignment methods or alignment-free phylogenetic approaches
- Studying the impact of alignment errors on tree inference
- Generating realistic benchmarks that mirror empirical data
- Comparing phylogenetic methods under challenging conditions
- Simulating protein evolution (indels are common in proteins)

**Skip indel simulation when:**
- Quick testing or prototyping (substitutions-only is faster)
- Studying tree inference methods assuming perfect alignment
- Computational resources are limited
- Alignment quality is not a concern for the research question

### Biological Rationale

Insertions and deletions are fundamental evolutionary processes occurring through:
- **Replication slippage**: Most common, produces short (1-2 bp/aa) indels
- **Recombination errors**: Produces medium-length indels
- **Transposon activity**: Produces longer insertions (rare)
- **Unequal crossing over**: Produces structural variation (very rare)

**Key empirical observations:**
1. **Deletions are more common than insertions**: ~3:1 ratio (deletion:insertion)
2. **Most indels are very short**: 60-80% are single-character indels
3. **Indel length follows power-law distribution**: Longer indels become exponentially rarer
4. **Context-dependent variation**: Indel rates vary across genomic regions

### Indel Parameters (Empirical)

The pipeline uses empirically-derived parameters from multiple studies (Benner et al. 1993; Cartwright 2009; Ly-Trong et al. 2022):

**Default Parameters:**
```bash
--indel 0.03,0.09       # Insertion rate: 0.03, Deletion rate: 0.09 (relative to substitution rate)
--indel-size POW 1.7 50 # Zipfian (power-law) distribution, exponent 1.7, max length 50
```

**Important Note on DNA vs Protein Sequences:**

The same indel parameters are applied to both DNA and protein sequences in this pipeline. This is a simplification appropriate for tree inference benchmarking:

- **For protein sequences**: The parameters (rate 0.12 total, Zipfian exponent 1.7) are empirically derived from protein evolution studies (Benner et al. 1993) and are biologically accurate.

- **For DNA sequences**: The total rate of 0.12 (12 indels per 100 substitutions) matches empirical rates for non-coding DNA. However, for coding DNA, indels would ideally be constrained to multiples of 3 nucleotides (whole codons) to avoid frameshifts. AliSim does not support this constraint, and implementing codon-aware indels would require a different simulator (e.g., INDELible with codon models).

- **Justification**: Since the primary use case is protein sequences and the focus is tree inference rather than realistic coding sequence evolution, using the same simplified indel model for both sequence types provides consistency and simplicity. The empirical rates are reasonable approximations for both contexts.

For studies specifically focused on modeling coding DNA evolution with frameshifts, consider using codon-aware simulators like INDELible.

**Parameter Explanation:**

1. **Indel Rates (--indel INS,DEL)**
   - **INS**: Insertion rate relative to substitution rate
   - **DEL**: Deletion rate relative to substitution rate
   - **Ratio**: DEL/INS ≈ 3.0 (deletions ~3× more common)
   - **Total indel rate**: 0.03 + 0.09 = 0.12 (12% relative to substitution rate)
   - **Biological interpretation**: For every 100 substitutions, expect ~3 insertions and ~9 deletions

2. **Indel Length Distribution (--indel-size)**

   AliSim supports multiple distributions:
   - **POW a max**: Zipfian/Power-law distribution (recommended)
     - `a`: Power-law exponent (1.7 = empirical estimate from protein data)
     - `max`: Maximum indel length (50 = reasonable biological upper bound)
   - **GEO p**: Geometric distribution
     - `p`: Probability parameter (0 < p ≤ 1)
   - **NB r p**: Negative Binomial distribution
     - `r`: Number of successes, `p`: probability
   - **LAV a max**: Lavalette distribution
     - Similar to Zipfian but with different tail behavior

**Zipfian Distribution (Power-Law):**

The Zipfian distribution with exponent 1.7 produces:
- ~60-70% single-character indels (1 bp or 1 aa)
- ~15-20% two-character indels
- ~5-10% three-character indels
- <5% indels of length ≥4
- Exponentially decreasing probability for longer indels
- Maximum length cap prevents unrealistically long indels

**Why exponent 1.7?**
- Estimated from empirical studies of protein evolution (Benner et al. 1993)
- Validated across pseudogenes and non-coding regions (Gu & Li 1995)
- Matches observed distributions in fungal genomes (Rasmussen & Kellis 2011)
- Standard in phylogenetic simulation studies (Cartwright 2009)

### Parameter Recommendations by Use Case

The pipeline provides convenient presets via the `--indel-model` option:

**Preset Models (Quick Selection):**

```bash
# No indels - substitutions only
--indel-model noindels
# Equivalent to: (indels disabled)

# Realistic (default) - empirically-derived parameters
--indel-model realistic
# Equivalent to: --indel 0.03,0.09 --indel-size POW{1.7/50}

# Conservative - moderate indel effects
--indel-model conservative
# Equivalent to: --indel 0.02,0.04 --indel-size POW{1.7/20}

# High rate - stress testing
--indel-model highrate
# Equivalent to: --indel 0.06,0.18 --indel-size POW{1.7/50}
```

**Manual Parameter Specification:**

For full control, specify rates and distribution directly:

**1. Default/Realistic (Recommended):**
```bash
--indel 0.03,0.09 --indel-size POW 1.7 50
```
Use for: Standard benchmarking, realistic simulations

**2. Conservative (Less impact):**
```bash
--indel 0.02,0.04 --indel-size POW 1.7 20
```
Use for: Testing with moderate indel effects, shorter sequences

**3. High Indel Rate (Challenging):**
```bash
--indel 0.06,0.18 --indel-size POW 1.7 50
```
Use for: Stress-testing methods, highly divergent sequences

**4. Short Indels Only:**
```bash
--indel 0.03,0.09 --indel-size POW 1.7 10
```
Use for: Focusing on micro-indels, reducing alignment complexity

**5. Geometric Distribution (Alternative):**
```bash
--indel 0.03,0.09 --indel-size GEO 0.7
```
Use for: Simpler model, faster decay of indel lengths

### Alignment Considerations

**Important**: When indels are simulated, sequences have different lengths and AliSim produces **true alignments** with gap characters ("-").

**Two Approaches:**

**Approach 1: Use True Alignment (Current Pipeline Default)**
- AliSim outputs sequences with gaps showing true evolutionary history
- PhyML directly accepts gapped alignments
- **Advantages**: Tests tree inference methods assuming perfect alignment
- **Use when**: Studying tree inference methods independently of alignment quality

**Approach 2: Re-align Sequences (Optional, Not Implemented)**
- Remove gaps from AliSim output to get unaligned sequences
- Run alignment software (MAFFT, MUSCLE, Clustal, etc.)
- Infer trees from inferred alignments
- **Advantages**: Tests full pipeline including alignment errors
- **Use when**: Studying how alignment errors affect tree inference
- **Implementation**: Would require adding alignment step before PhyML

**Current Pipeline Choice:**

The pipeline uses **Approach 1** (true alignments) because:
1. Focuses on tree inference methods, not alignment methods
2. Faster execution (no alignment step needed)
3. Clean separation of alignment quality from tree inference quality
4. Standard practice in tree inference benchmarking (Molloy & Warnow 2020)

**Unaligned FASTA Files:**

AliSim automatically generates **`alignment.unaligned.fa`** files containing unaligned sequences with all gap characters removed. This file provides:
- **FASTA format**: Sequences with `>` headers for easy processing
- **No gaps**: All `-` and space characters removed from sequences
- **Unaligned sequences**: Ready for testing alignment software (MAFFT, MUSCLE, etc.)

**Use Cases:**
1. **Alignment benchmarking**: Use as input to alignment software and compare results to `alignment_TRUE.phy`
2. **Alignment-free methods**: Test phylogenetic methods that work directly on unaligned sequences
3. **Full pipeline testing**: Run alignment → tree inference → compare to true tree
4. **Format conversion**: FASTA is more widely supported than PHYLIP for many tools

**Example workflow for testing alignment methods:**
```bash
# 1. Run pipeline to generate sequences
./simulate_pipeline.sh -n 10 -l 12 --indel-model realistic

# 2. Align the unaligned sequences (example with MAFFT)
mafft alignment.unaligned.fa > inferred_alignment.fasta

# 3. Compare inferred alignment to true alignment
# (convert alignment_TRUE.phy to FASTA first if needed)

# 4. Infer tree from inferred alignment
phyml -i inferred_alignment.fasta ...

# 5. Compare inferred tree to ml_gene_tree.nwk
```

**To test alignment methods without using the nogaps file**, you would need to:
1. Extract unaligned sequences from AliSim output (remove gaps)
2. Run alignment software (e.g., MAFFT, MUSCLE)
3. Compare inferred alignment to true alignment
4. Infer trees from inferred alignment

### Literature Support

**Indel Model Development:**
- **Thorne-Kishino-Felsenstein (TKF91)**: Classic single-character indel model (Thorne et al. 1991)
- **TKF92**: Extended model allowing multi-character indels (Thorne et al. 1992)
- **Dawg**: Early simulator with power-law indel lengths (Cartwright 2005, 2009)
- **INDELible**: Flexible simulator with multiple indel distributions (Fletcher & Yang 2009)
- **AliSim**: Modern, fast simulator with empirical parameterization (Ly-Trong et al. 2022)

**Empirical Studies:**
- **Benner et al. (1993)**: Power-law distribution in proteins, exponent ~1.7
- **Gu & Li (1995)**: Indel rates in pseudogenes and coding sequences
- **Cartwright (2009)**: Zipfian distribution in sequence evolution
- **Rasmussen & Kellis (2011)**: Fungal genome indel rates (used in Molloy & Warnow 2020)

**Benchmarking Studies:**
- **Molloy & Warnow (2020)**: TreeShrink benchmarking (basis for our substitution parameters)
- **Zhu et al. (2025)**: Deep learning phylogenetic inference with indels
- **Fletcher & Yang (2009)**: INDELible validation across multiple distributions

### Example Usage

**Using Preset Models (Recommended):**
```bash
# No indels - substitutions only (fastest)
./simulate_pipeline.sh -n 10 -l 12 --indel-model noindels

# Realistic indels (empirically-derived parameters)
./simulate_pipeline.sh -n 10 -l 12 --indel-model realistic

# Conservative indels (lower rates, shorter lengths)
./simulate_pipeline.sh -n 10 -l 12 --indel-model conservative

# High indel rate (stress testing)
./simulate_pipeline.sh -n 10 -l 12 --indel-model highrate
```

**Manual Parameter Specification:**
```bash
# Enable indels with default empirical parameters
./simulate_pipeline.sh -n 10 -l 12 --indel 0.03,0.09

# Custom indel parameters
./simulate_pipeline.sh -n 10 -l 12 --indel 0.02,0.04 --indel-size "POW 1.7 20"

# High indel rate for stress testing
./simulate_pipeline.sh -n 10 -l 12 --indel 0.06,0.18 --indel-size "POW 1.7 50"

# Geometric distribution
./simulate_pipeline.sh -n 10 -l 12 --indel 0.03,0.09 --indel-size "GEO 0.7"
```

### Output with Indels

When indels are enabled:
- **Sequence lengths vary** due to insertions and deletions
- **Alignments contain gaps** ("-" characters) showing true evolutionary history
- **PhyML handles gapped alignments** directly without re-alignment
- **ML trees estimated** from true alignments including gap patterns
- All other pipeline behavior remains unchanged

### Validation

The indel model can be validated by:
1. **Indel count distributions**: Should match power-law with exponent 1.7
2. **Indel length distributions**: Should match specified distribution
3. **Insertion vs deletion ratio**: Should match specified rates (default 1:3)
4. **Sequence length variation**: Should reflect cumulative indel effects
5. **Gap patterns in alignment**: Should show realistic clustering

## Output Structure

For each replicate and configuration:

```
simulation_output/
├── simulation_summary.txt                # Complete parameter summary and file statistics
├── species_trees/
│   └── lf12_replicate1.nex               # Species tree for Phase 2 (Nexus format)
├── gene_trees/
│   └── leaves12_dl1e-10_ps1e7/
│       └── replicate_1/
│           └── 1/
│               ├── s_tree.trees          # Species tree (from SimPhy)
│               ├── l_trees.trees         # Locus tree (after DL)
│               └── g_trees1.trees        # Gene tree (after ILS, branches in subs/site)
├── dna/
│   └── leaves12_dl1e-10_ps1e7/
│       └── replicate_1/
│           ├── gene_tree.nwk             # Gene tree from Phase 2 (branches in subs/site)
│           ├── alignment_TRUE.phy        # True alignment (from AliSim)
│           ├── alignment.unaligned.fa    # Unaligned sequences (auto-generated by AliSim)
│           ├── alignment_TRUE.phy_phyml_tree.txt  # ML tree
│           ├── alignment_TRUE.phy_phyml_stats.txt # PhyML statistics
│           └── ml_gene_tree.nwk          # Copied ML tree
└── protein/
    └── leaves12_dl1e-10_ps1e7/
        └── replicate_1/
            ├── gene_tree.nwk             # Gene tree from Phase 2 (branches in subs/site)
            ├── alignment_TRUE.phy        # True alignment (from AliSim)
            ├── alignment.unaligned.fa    # Unaligned sequences (auto-generated by AliSim)
            ├── alignment_TRUE.phy_phyml_tree.txt  # ML tree
            ├── alignment_TRUE.phy_phyml_stats.txt # PhyML statistics
            └── ml_gene_tree.nwk          # Copied ML tree
```

### Simulation Summary File

The pipeline automatically generates `simulation_summary.txt` in the output directory containing:

**Configuration Documentation:**
- All command-line parameters used
- Species tree, gene tree, and sequence parameters
- Indel simulation settings (if enabled)
- Gene tree filtering settings (if enabled)
- Timestamp, hostname, and working directory

**File Statistics:**
- Count of species trees generated
- Count of gene trees generated
- Count of DNA alignments and ML trees
- Count of protein alignments and ML trees

**Additional Information:**
- Software versions (SimPhy, IQ-TREE 2, PhyML)
- Complete directory structure
- Configuration matrix listing all parameter combinations

This file provides complete documentation for reproducibility and sharing simulation details with collaborators.

## Simulation Scale

With default parameters:
- **Leaf-set sizes**: 2 (12 and 20 leaves)
- **DL rates**: 3 levels
- **Population sizes**: 2 levels
- **Replicates per configuration**: 100
- **Total configurations**: 2 × 3 × 2 = 12
- **Total datasets per sequence type**: 1,200
- **Total datasets (DNA + Protein)**: 2,400

## Pipeline Improvements (2025-12-23)

The pipeline was significantly improved to use SimPhy end-to-end:

### Previous Approach (TreeSim + SimPhy)
- Phase 1: TreeSim generated species trees
- Required scaling branch lengths (years → millions of years)
- Required scaling all rates by 10⁶
- Required regex manipulation to remove root branches
- Prone to scaling errors and tree format issues

### Current Approach (SimPhy Only)
- Phase 1: **SimPhy generates species trees directly**
- Phase 2: SimPhy generates gene trees (same tool, consistent behavior)
- No scaling needed - all parameters in years
- SimPhy outputs gene tree branch lengths directly in **substitutions per site**
- No tree manipulation required
- Simpler, more reliable, fewer sources of error

### Bug Fixes Applied
1. **Path mismatch**: Phase 2 and Phase 3 used inconsistent replicate directory names
2. **Incorrect scaling removal**: Removed obsolete branch length scaling (×0.0004) from Phase 3/4
3. **Octal number errors**: Fixed arithmetic with zero-padded numbers (01-09)

**⚠️ WARNING**: Any simulations generated with the old TreeSim-based pipeline should be regenerated.

See `BUGFIX_DL_AND_SEED.md` for documentation of earlier fixes.

## References

### Core Pipeline

- **Molloy, E.K. and Warnow, T. (2020).** TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees. BMC Genomics, 21(Suppl 2):332.
- **Rasmussen, M.D. and Kellis, M. (2011).** A Bayesian approach for fast and accurate gene tree reconstruction. Molecular Biology and Evolution, 28(1):273-290.
- **Mallo, D., De Oliveira Martins, L., and Posada, D. (2016).** SimPhy: Phylogenomic simulation of gene, locus, and species trees. Systematic Biology, 65(2):334-344.
- **Ly-Trong, N. et al. (2022).** AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era. Molecular Biology and Evolution, 39(5):msac092.
- **Minh, B.Q. et al. (2020).** IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5):1530-1534.
- **Guindon, S. et al. (2010).** New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology, 59(3):307-321.

### Indel Modeling

**Model Development:**
- **Thorne, J.L., Kishino, H., and Felsenstein, J. (1991).** An evolutionary model for maximum likelihood alignment of DNA sequences. Journal of Molecular Evolution, 33(2):114-124.
- **Thorne, J.L., Kishino, H., and Felsenstein, J. (1992).** Inching toward reality: An improved likelihood model of sequence evolution. Journal of Molecular Evolution, 34(1):3-16.
- **Cartwright, R.A. (2005).** DNA assembly with gaps (Dawg): simulating sequence evolution. Bioinformatics, 21(Suppl 3):iii31-iii38.
- **Cartwright, R.A. (2009).** Problems and solutions for estimating indel rates and length distributions. Molecular Biology and Evolution, 26(2):473-480.
- **Fletcher, W. and Yang, Z. (2009).** INDELible: A flexible simulator of biological sequence evolution. Molecular Biology and Evolution, 26(8):1879-1888.

**Empirical Studies:**
- **Benner, S.A., Cohen, M.A., and Gonnet, G.H. (1993).** Empirical and structural models for insertions and deletions in the divergent evolution of proteins. Journal of Molecular Biology, 229(4):1065-1082.
- **Gu, X. and Li, W.H. (1995).** The size distribution of insertions and deletions in human and rodent pseudogenes suggests the logarithmic gap penalty for sequence alignment. Journal of Molecular Evolution, 40(4):464-473.

**Recent Benchmarking:**
- **Zhu, T. et al. (2025).** A critical evaluation of deep-learning based phylogenetic inference. Journal of Genetics and Genomics, 52(1):in press.

## Command-Line Usage

The pipeline can be run with customizable parameters:

```bash
# Default: 100 replicates, 12 and 20 leaves, both DNA and protein
./simulate_pipeline.sh

# Custom configuration
./simulate_pipeline.sh -n 50 -l 12,20,50 -t both -o my_output -s 42

# Two-stage workflow: alignments first, ML trees later
# Stage 1: Generate alignments only (skip PhyML, ~3x faster)
./simulate_pipeline.sh -n 100 -l 12 -o my_output -m

# Stage 2: Infer ML trees from existing alignments
./simulate_pipeline.sh -o my_output -n 100 -l 12 -i

# Options:
#   -n NUM               Number of replicates per configuration (default: 100)
#   -l LEAVES            Comma-separated leaf counts (default: 12,20)
#   -t TYPE              Sequence type: dna, protein, or both (default: both)
#   -o DIR               Output directory (default: simulation_output)
#   -s SEED              Random seed (default: 22)
#   -r MAX               Retry if duplicate sequences found, up to MAX attempts (default: 0)
#   -f NUM               Filter gene trees: require exactly NUM leaves (default: 0 = no filtering)
#   -m                   Skip ML tree estimation (generate alignments only)
#   -i                   Infer-only mode: run only PhyML on existing alignments
#   --dna-length NUM     DNA alignment length in base pairs (default: 1000)
#   --protein-length NUM Protein alignment length in amino acids (default: 333)
#   --indel-model MODEL  Indel preset: noindels, realistic, conservative, or highrate
#   --indel INS,DEL      Enable indel simulation with specified rates (e.g., 0.03,0.09)
#   --indel-size DIST    Indel length distribution (default: POW 1.7 50)
```

### Two-Stage Workflow

The pipeline supports a two-stage workflow for flexibility and efficiency:

**Stage 1: Simulation and Alignment Generation** (`-m` flag)
- Runs Phases 1-4 but **skips PhyML** tree inference
- Generates species trees, gene trees, and sequence alignments
- Approximately **3× faster** than full pipeline
- Useful for quickly generating alignment datasets

**Stage 2: ML Tree Inference** (`-i` flag)
- **Skips** Phases 1-2 (species and gene tree simulation)
- **Skips** sequence generation in Phases 3-4
- **Runs only PhyML** on existing `alignment_TRUE.phy` files
- Generates ML trees and PhyML statistics
- Must use **same parameters** (`-n`, `-l`) as Stage 1

**Use Cases:**
1. **Fast testing**: Generate alignments quickly, infer trees only if needed
2. **Batch processing**: Generate alignments once, distribute ML inference across compute nodes
3. **Parameter exploration**: Generate alignments once, experiment with different ML inference settings

### Gene Tree Filtering by Leaf Count

Due to duplication and loss events in Phase 2, gene trees can have varying numbers of leaves. The `-f NUM` option filters gene trees to ensure all replicates have exactly NUM leaves.

**Example:**
```bash
# Require all gene trees to have exactly 10 leaves
./simulate_pipeline.sh -n 100 -l 12 -f 10
```

**Filtering Process:**
1. SimPhy generates a gene tree with duplication/loss and ILS
2. The pipeline counts the number of leaves in the gene tree
3. If the tree has exactly NUM leaves, it's accepted for Phases 3-4
4. If not, SimPhy regenerates with a different seed (up to 1000 attempts)
5. Progress messages indicate how many retries were needed

**Technical Details:**
- Leaf counting uses a simple algorithm: count commas in Newick format + 1
- Each retry uses a different random seed: `seed + retry_count × 1,000,000`
- Maximum 1000 retries per replicate (configurable via `MAX_GENE_TREE_RETRIES`)
- If maximum retries reached, the last tree is accepted with a warning

**Example output:**
```
Configuration: leaves12_dl1e-10_ps1e7
  Replicate 1: accepted after 15 retries (10 leaves)
  Replicate 2: accepted after 3 retries (10 leaves)
  ...
```

**Use Cases:**
- **Consistent tree sizes**: Ensure all datasets have the same number of taxa
- **Fair comparison**: Remove variability due to different tree sizes
- **Controlled experiments**: Study effects of ILS and DL while maintaining constant tree size
- **Filtering artifacts**: Exclude gene trees with extreme duplication/loss outcomes

### Custom Alignment Lengths

The pipeline allows customization of alignment lengths for both DNA and protein sequences via command-line options.

**Example:**
```bash
# Custom protein alignment length
./simulate_pipeline.sh -n 100 -l 12 --protein-length 500

# Custom DNA alignment length
./simulate_pipeline.sh -n 100 -l 12 --dna-length 2000

# Both custom lengths
./simulate_pipeline.sh -n 100 -l 12 --dna-length 1500 --protein-length 500

# Longer alignments with indels
./simulate_pipeline.sh -n 100 -l 12 --protein-length 500 --indel-model realistic
```

**Technical Details:**
- DNA and protein alignment lengths are set independently
- The `--length` parameter is passed directly to AliSim
- **Without indels**: Alignment length equals the specified value exactly
- **With indels**: Final alignment length varies due to insertions and deletions
  - Expected length ≈ specified length × (1 + insertion_rate - deletion_rate)
  - Example: 1000 bp with `--indel 0.03,0.09` → ~940 bp on average

**Use Cases:**
- **Longer alignments**: Improve phylogenetic signal for difficult inference problems
- **Shorter alignments**: Faster simulations for testing and prototyping
- **Variable lengths**: Test method robustness across different data sizes
- **Protein-specific lengths**: Simulate realistic protein domain lengths (e.g., 200-500 aa)

See `README.md` for complete usage documentation.
