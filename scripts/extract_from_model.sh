#!/usr/bin/env bash
set -euo pipefail

# --- Required positional args ---
DATASET="$1"
MODEL="$2"
INPUT_FASTA="$3"
OUTPUT_CSV="$4"
shift 4  # shift away required args, remaining are optional kwargs

# --- Optional keyword args ---
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-size|--layer|--max-length|--precision|--pooling|--pca-dims)
            EXTRA_ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            echo "[WARN] Unknown argument: $1" >&2
            shift
            ;;
    esac
done

# --- GPU slot management ---
SLOT=$(python scripts/acquire_slot.py)
echo "[GPU slot] $SLOT"
export CUDA_VISIBLE_DEVICES="$SLOT"

cleanup() {
    python scripts/release_slot.py "$SLOT"
}
trap cleanup EXIT

# --- Output directory (derived safely) ---
OUTDIR=$(dirname "$OUTPUT_CSV")
mkdir -p "$OUTDIR"

# --- Run extractor ---
python "scripts/extract_${MODEL}.py" \
    --input "$INPUT_FASTA" \
    --output "$OUTPUT_CSV" \
    --gpu-id 0 \
    "${EXTRA_ARGS[@]}"