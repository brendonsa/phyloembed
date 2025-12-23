#!/usr/bin/env python3
import argparse
import os
import subprocess
import tempfile
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", required=True, help="results/<dataset> directory")
    parser.add_argument("--rscript", required=True, help="Path to compare_trees.R")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument(
        "--ref-name",
        default="tree.nwk",
        help="Reference tree filename inside dataset-dir (default: tree.nwk)",
    )
    args = parser.parse_args()

    dataset_dir = args.dataset_dir
    ref_name = args.ref_name
    ref_path = os.path.join(dataset_dir, ref_name)

    if not os.path.exists(ref_path):
        raise SystemExit(f"Reference tree not found: {ref_path}")

    # non-recursive: only .nwk directly in this directory
    nwk_files = [
        os.path.join(dataset_dir, f)
        for f in os.listdir(dataset_dir)
        if f.endswith(".nwk")
    ]

    others = [p for p in nwk_files if os.path.basename(p) != ref_name]

    rows = []
    for tree_path in others:
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            tmp_path = tmp.name

        subprocess.run(
            [
                "Rscript",
                args.rscript,
                "--t1",
                ref_path,
                "--t2",
                tree_path,
                "--out",
                tmp_path,
            ],
            check=True,
        )

        df = pd.read_csv(tmp_path)  # columns: metric, value
        os.remove(tmp_path)

        # pivot to wide: one row, columns = metrics
        metrics = df.set_index("metric")["value"].to_dict()
        metrics["tree"] = os.path.basename(tree_path)
        rows.append(metrics)

    if rows:
        out_df = pd.DataFrame(rows)
        # nice column order if metrics are present
        cols = ["tree", "RF_raw", "RF_normalized", "Weighted_RF", "BranchScore_KF", "JRF", "JRF_normalized", "Nye", "Nye_normalized","SMI","SMI_normalized"]
        out_df = out_df[[c for c in cols if c in out_df.columns]]
    else:
        out_df = pd.DataFrame(columns=["tree","RF_raw", "RF_normalized", "Weighted_RF", "BranchScore_KF", "JRF", "JRF_normalized", "Nye", "Nye_normalized","SMI","SMI_normalized"])

    out_df.sort_values(by='RF_normalized', inplace=True)

    out_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
