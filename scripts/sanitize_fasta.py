#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--map", default=None)
    args = ap.parse_args()

    in_fa = Path(args.input)
    out_fa = Path(args.output)
    map_file = Path(args.map) if args.map else out_fa.with_suffix(".labels.tsv")

    labels = []
    seqs = []

    with open(in_fa, "r", encoding="utf-8", errors="strict") as f:
        cur = None
        for line in f:
            if line.startswith(">"):
                cur = line[1:].strip()
                labels.append(cur)
                seqs.append([])
            else:
                seqs[-1].append(line.rstrip())

    ids = [f"t{i:05d}" for i in range(1, len(labels) + 1)]

    with open(map_file, "w", encoding="utf-8") as m:
        for old, new in zip(labels, ids):
            m.write(f"{old}\t{new}\n")

    with open(out_fa, "w", encoding="utf-8") as out:
        for new, s in zip(ids, seqs):
            out.write(f">{new}\n")
            out.write("\n".join(s) + "\n")

if __name__ == "__main__":
    main()
