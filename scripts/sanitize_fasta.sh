#!/usr/bin/env bash
set -euo pipefail

input="$1"
output="$2"

mkdir -p "$(dirname "$output")"

awk '
  /^>/{
    h = substr($0, 2)
    gsub(/[ =]/, "_", h)          # replace space and '=' with '_'
    if (length(h) > 64) h = substr(h, 1, 64)
    print ">" h
    next
  }
  { print }
' "$input" > "$output"
