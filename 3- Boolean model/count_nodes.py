#!/usr/bin/env python3
"""
count_bnet_folder.py

Applies the total / constant / non-constant node counting to ALL .bnet
files in a folder (e.g. your 1000 models), and outputs a per-file summary
+ global statistics (min/max/mean).

Usage:
    python count_bnet_folder.py /path/to/models_folder
    python count_bnet_folder.py /path/to/models_folder --csv results.csv
    python count_bnet_folder.py /path/to/models_folder --ext .bnet
"""
import argparse
import csv
import os
import sys


def count_bnet(model_file):
    n_non_constant = 0
    n_constant = 0
    with open(model_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("targets"):
                continue
            if "," not in line:
                continue
            node, rule = [x.strip() for x in line.split(",", 1)]
            if rule in {"0", "1"}:
                n_constant += 1
            else:
                n_non_constant += 1
    total = n_non_constant + n_constant
    return total, n_non_constant, n_constant


def main():
    parser = argparse.ArgumentParser(
        description="Count total / constant / non-constant nodes for all .bnet files in a folder"
    )
    parser.add_argument("folder", help="Folder containing the .bnet files (e.g. the 1000 models)")
    parser.add_argument("--ext", default=".bnet", help="File extension to analyze (default: .bnet)")
    parser.add_argument("--csv", default=None, help="Path to export the results as CSV")
    args = parser.parse_args()

    if not os.path.isdir(args.folder):
        print(f"Error: '{args.folder}' is not a valid folder.", file=sys.stderr)
        sys.exit(1)

    files = sorted(
        f for f in os.listdir(args.folder) if f.endswith(args.ext)
    )

    if not files:
        print(f"No '{args.ext}' files found in {args.folder}", file=sys.stderr)
        sys.exit(1)

    results = []
    for fname in files:
        fpath = os.path.join(args.folder, fname)
        total, n_non_const, n_const = count_bnet(fpath)
        results.append(
            {"file": fname, "total": total, "n_non_constant": n_non_const, "n_constant": n_const}
        )

    # Per-file display
    print(f"{'File':<30} {'Total':>8} {'Non-const':>10} {'Const':>8}")
    print("-" * 60)
    for r in results:
        print(f"{r['file']:<30} {r['total']:>8} {r['n_non_constant']:>10} {r['n_constant']:>8}")

    # Global statistics
    totals = [r["total"] for r in results]
    non_consts = [r["n_non_constant"] for r in results]
    consts = [r["n_constant"] for r in results]

    print("-" * 60)
    print(f"Number of models analyzed : {len(results)}")
    print(
        f"Total nodes        -> min={min(totals)} max={max(totals)} "
        f"mean={sum(totals)/len(totals):.2f}"
    )
    print(
        f"Non-constant nodes -> min={min(non_consts)} max={max(non_consts)} "
        f"mean={sum(non_consts)/len(non_consts):.2f}"
    )
    print(
        f"Constant nodes     -> min={min(consts)} max={max(consts)} "
        f"mean={sum(consts)/len(consts):.2f}"
    )

    # If all models have the same total number of nodes (common case
    # within a single BoNesis inference run), display it clearly
    if len(set(totals)) == 1:
        print(f"\n=> All models have the same total number of nodes: {totals[0]}")
    if len(set(non_consts)) == 1:
        print(f"=> All models have the same number of non-constant nodes: {non_consts[0]}")

    # CSV export
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["file", "total", "n_non_constant", "n_constant"])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSummary exported to: {args.csv}")


if __name__ == "__main__":
    main()
