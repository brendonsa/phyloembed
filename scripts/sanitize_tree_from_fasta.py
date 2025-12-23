#!/usr/bin/env python3
import argparse

def load_map(path):
    m = {}
    with open(path) as f:
        for line in f:
            a, b = line.rstrip().split("\t")
            m[a] = b
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    m = load_map(args.map)

    with open(args.tree) as f:
        s = f.read()

    for k, v in m.items():
        s = s.replace(k, v)

    with open(args.output, "w") as out:
        out.write(s)

if __name__ == "__main__":
    main()
