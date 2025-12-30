# Quick Start Guide

## Prerequisites

Install required software:

```bash
# Download and install SimPhy
wget https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz
tar -xzf SimPhy_1.0.2.tar.gz
chmod +x SimPhy_1.0.2/bin/simphy_lnx64

# The pipeline auto-detects SimPhy in SimPhy_1.0.2/bin/
# Or optionally add to PATH:
# export PATH=$PATH:$(pwd)/SimPhy_1.0.2/bin

# Install IQ-TREE 2 (includes AliSim for sequence simulation)
# IMPORTANT: Version 2.2.0+ required for AliSim
# Repository versions may be outdated - manual install recommended
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz
tar -xzf iqtree-2.3.6-Linux-intel.tar.gz
sudo cp iqtree-2.3.6-Linux-intel/bin/iqtree2 /usr/local/bin/

# Install PhyML (for ML tree estimation)
# Ubuntu/Debian:
sudo apt-get install phyml

# Verify installations
simphy_lnx64 -h 2>&1 | head -10
iqtree2 --version
phyml --version
```

## Running the Pipeline

### Option 1: Bash Script (Recommended)

```bash
# Default: 100 replicates, 12 and 20 leaves, both DNA and protein
./simulate_pipeline.sh

# Quick test: 10 replicates, 12 leaves only, both types
./simulate_pipeline.sh -n 10 -l 12

# Run DNA only (faster)
./simulate_pipeline.sh -t dna

# Run protein only
./simulate_pipeline.sh -t protein

# Custom configuration
./simulate_pipeline.sh -n 50 -l 12,20,50 -t both -o my_output -s 42

# Show all options
./simulate_pipeline.sh -h
```

### Option 2: Python Script

```bash
# Default configuration
./simulate_pipeline.py

# Quick test
./simulate_pipeline.py -n 10 -l 12

# DNA only
./simulate_pipeline.py --type dna

# Protein only with custom leaves
./simulate_pipeline.py -t protein -l 12 20 50

# Full custom configuration
./simulate_pipeline.py -n 50 -l 12 20 -t both -o my_output -s 42

# Show all options
./simulate_pipeline.py -h
```

### Command-Line Options

Both scripts support the same options:

| Option | Bash | Python | Description | Default |
|--------|------|--------|-------------|---------|
| Replicates | `-n NUM` | `-n, --replicates NUM` | Number of replicates | 100 |
| Leaves | `-l LEAVES` | `-l, --leaves LEAVES` | Leaf counts | 12,20 (bash) or 12 20 (python) |
| Type | `-t TYPE` | `-t, --type TYPE` | dna, protein, or both | both |
| Output | `-o DIR` | `-o, --output DIR` | Output directory | simulation_output |
| Seed | `-s SEED` | `-s, --seed SEED` | Random seed | 22 |
| Filter | `-f NUM` | N/A | Require gene trees with exactly NUM leaves | 0 (no filter) |
| Skip ML | `-m` | N/A | Skip PhyML (alignments only) | - |
| Infer-only | `-i` | N/A | Run only PhyML on existing alignments | - |
| Help | `-h` | `-h, --help` | Show help | - |

## What Gets Generated

For **each configuration** (12 total by default):
- 100 species trees
- 100 gene trees (with duplication/loss and ILS)
- 100 DNA MSAs + ML trees
- 100 protein MSAs + ML trees

**Total:** 1,200 datasets for DNA + 1,200 datasets for proteins

## Output Structure

```
simulation_output/
├── species_trees/       # Species trees (Newick and Nexus)
├── gene_trees/          # True gene trees from SimPhy
├── dna/                 # DNA alignments and ML trees
└── protein/             # Protein alignments and ML trees
```

## Configuration Matrix

| Leaves | DL Rate | Pop Size | ILS Level |
|--------|---------|----------|-----------|
| 12     | 1e-10   | 1e7      | Low       |
| 12     | 1e-10   | 5e7      | Medium    |
| 12     | 2e-10   | 1e7      | Low       |
| 12     | 2e-10   | 5e7      | Medium    |
| 12     | 5e-10   | 1e7      | Low       |
| 12     | 5e-10   | 5e7      | Medium    |
| 20     | 1e-10   | 1e7      | Low       |
| 20     | 1e-10   | 5e7      | Medium    |
| 20     | 2e-10   | 1e7      | Low       |
| 20     | 2e-10   | 5e7      | Medium    |
| 20     | 5e-10   | 1e7      | Low       |
| 20     | 5e-10   | 5e7      | Medium    |

## Estimated Runtime

- **Phase 1** (Species trees): ~5-10 minutes
- **Phase 2** (Gene trees): ~30-60 minutes
- **Phase 3** (DNA sequences): ~2-4 hours
- **Phase 4** (Protein sequences): ~2-4 hours

**Total:** ~5-9 hours for complete pipeline with default parameters

## Disk Space

Expect to use approximately **5-10 GB** for the complete output with default parameters.

## Customization Examples

### Smaller Test Run

```bash
# Quick test - only 10 datasets
./simulate_pipeline.sh -n 10 -l 12

# Or with Python
./simulate_pipeline.py -n 10 -l 12
```

This generates only 10 datasets (2 configs: low/medium ILS, but only 12 leaves and default DL rates).

### DNA Only (Faster)

```bash
# Generate DNA sequences only
./simulate_pipeline.sh -t dna

# Or with Python
./simulate_pipeline.py --type dna
```

### Protein Only

```bash
# Generate protein sequences only
./simulate_pipeline.sh -t protein

# Or with Python
./simulate_pipeline.py --type protein
```

### Large Phylogeny (50 taxa)

```bash
# Simulate trees with 50 leaves
./simulate_pipeline.sh -l 50 -n 50

# Or with Python
./simulate_pipeline.py -l 50 -n 50
```

### Multiple Leaf Sizes

```bash
# Bash: comma-separated
./simulate_pipeline.sh -l 12,20,50

# Python: space-separated
./simulate_pipeline.py -l 12 20 50
```

### Two-Stage Workflow (Alignments First, ML Trees Later)

For faster testing, you can generate alignments first, then infer ML trees later:

**Stage 1: Generate alignments only**
```bash
# Skip PhyML tree estimation (~3x faster)
./simulate_pipeline.sh -n 10 -l 12 -o my_test -m
```

This creates alignments but skips ML tree inference.

**Stage 2: Infer ML trees from existing alignments**
```bash
# Later: run only PhyML on existing alignments
./simulate_pipeline.sh -o my_test -n 10 -l 12 -i
```

This runs only PhyML on the existing `alignment_TRUE.phy` files, creating ML trees.

**Important:** Use the same `-n` and `-l` values in both stages!

### Gene Tree Filtering

Filter gene trees to have exactly a specific number of leaves (useful for controlling duplication/loss effects):

```bash
# Require all gene trees to have exactly 10 leaves
./simulate_pipeline.sh -n 100 -l 12 -f 10
```

This retries gene tree generation until a tree with exactly 10 leaves is found. Progress messages show the number of retries needed.

### Different Alignment Lengths (Advanced)

Edit the script directly:
```bash
ALIGNMENT_LENGTH_DNA=500
ALIGNMENT_LENGTH_PROTEIN=166
```

## Troubleshooting

### Check if software is installed

```bash
which simphy_lnx64
which iqtree2
which phyml
```

All should return paths. If not, install the missing software.

### Test installations

```bash
# SimPhy
simphy_lnx64 -h 2>&1 | head -5

# IQ-TREE 2 (should be 2.2.0+)
iqtree2 --version

# PhyML
phyml --version
```

### Permission denied

```bash
chmod +x simulate_pipeline.sh
chmod +x simulate_pipeline.py
```

## Next Steps

After running the pipeline:
1. Check `simulation_output/` directory
2. Explore the generated trees and alignments
3. Use the data for your phylogenetic analyses

## Getting Help

See `README.md` for detailed documentation and parameter explanations.
