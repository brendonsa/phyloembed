#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   iqtree_support.sh BASE_TREE BOOTSTRAP_TREES... OUTPUT_TREE
#
# Example:
#   iqtree_support.sh best.tree boots/*.tre best_with_bs.tree

if [ "$#" < 3 ]; then
    echo "Usage: $0 BASE_TREE BOOTSTRAP_TREES... OUTPUT_TREE" >&2
    exit 1
fi

BASE_TREE="$1"
shift

# Last arg = output tree
OUTPUT_TREE="${@: -1}"

# All remaining args (except last) = bootstrap tree files
BOOT_FILES=("${@:1:$#-1}")

# Concatenate all bootstrap trees into one file
BOOTCAT="${OUTPUT_TREE}.combined_boots.tre"
> "$BOOTCAT"
for f in "${BOOT_FILES[@]}"; do
    cat "$f" >> "$BOOTCAT"
done

PREFIX="${OUTPUT_TREE}.tmp"

# IQ-TREE: map support from bootstrap trees (-t) onto target tree (-sup)
iqtree2 \
    -t "$BOOTCAT" \
    -sup "$BASE_TREE" \
    -pre "$PREFIX" \
    -nt 1

# IQ-TREE writes: PREFIX.suptree
mv "${PREFIX}.suptree" "$OUTPUT_TREE"

# Clean up
rm -f "${PREFIX}".* "$BOOTCAT"
