#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   iqtree_consensus.sh BOOTSTRAP_TREES... OUTPUT_TREE
#
# Example:
#   iqtree_consensus.sh boots/*.tre consensus.tree

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 BOOTSTRAP_TREES... OUTPUT_TREE" >&2
    exit 1
fi

# last arg = output consensus tree
OUTPUT_TREE="${@: -1}"
# all previous args = bootstrap tree files
BOOT_FILES=("${@:1:$#-1}")

# combine all bootstrap trees into one file
BOOTCAT="${OUTPUT_TREE}.combined_boots.tre"
> "$BOOTCAT"
for f in "${BOOT_FILES[@]}"; do
    cat "$f" >> "$BOOTCAT"
done

PREFIX="${OUTPUT_TREE}.tmp"

# IQ-TREE 2.0.7 consensus:
# -t FILE   : set of input trees
# -con      : compute consensus tree -> .contree
# -minsup 0.5 : majority-rule consensus (0 = extended MR)
iqtree2 \
    -t "$BOOTCAT" \
    -con \
    -minsup 0.5 \
    -pre "$PREFIX" \
    -nt 1

# IQ-TREE writes PREFIX.contree
mv "${PREFIX}.contree" "$OUTPUT_TREE"

# clean up
rm -f "${PREFIX}".* "$BOOTCAT"
