#!/usr/bin/env python3
"""
count_attractors_folder.py - Count the attractors of ALL Boolean network
models (.bnet, BoolNet format) in a folder (e.g. your 1000 models),
and produce a summary table (one row per model).

Usage:
    python count_attractors_folder.py <folder>
    python count_attractors_folder.py <folder> --fixed-points
    python count_attractors_folder.py <folder> --from "a=0,b=1,c=0"
    python count_attractors_folder.py <folder> --csv summary.csv
    python count_attractors_folder.py <folder> --ext .bnet

Options:
    --fixed-points      Count fixed points only (instead of all attractors)
    --from STATE        Restrict to attractors reachable from STATE
                        Format: "node=value,node=value,..."
    --csv FILE          Write the summary table (model, count) to FILE
    --ext EXT           File extension to look for (default: .bnet)
"""
import sys
import os
import argparse
import pandas as pd


def parse_state(state_str):
    """Parse 'a=0,b=1,c=1' -> {'a': 0, 'b': 1, 'c': 1}"""
    result = {}
    for part in state_str.split(","):
        part = part.strip()
        if "=" not in part:
            print(f"Error: invalid format '{part}', expected 'node=value'", file=sys.stderr)
            sys.exit(1)
        key, val = part.split("=", 1)
        result[key.strip()] = int(val.strip())
    return result


def count_for_model(mpbn, bnet_file, fixed_points=False, reachable_from=None):
    """Load one .bnet file and return its attractor (or fixed-point) count."""
    mbn = mpbn.MPBooleanNetwork(bnet_file)

    kwargs = {}
    if reachable_from:
        kwargs["reachable_from"] = reachable_from

    if fixed_points:
        count = mbn.count_fixedpoints(**kwargs)
    else:
        count = mbn.count_attractors(**kwargs)

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Count attractors for all .bnet models in a folder using mpbn"
    )
    parser.add_argument("folder", help="Folder containing the .bnet model files (e.g. the 1000 models)")
    parser.add_argument(
        "--fixed-points", action="store_true",
        help="Count fixed points only (instead of all attractors)"
    )
    parser.add_argument(
        "--from", dest="init_state", default=None, metavar="STATE",
        help="Initial state, e.g.: 'a=0,b=1,c=0'"
    )
    parser.add_argument(
        "--ext", default=".bnet",
        help="File extension to look for (default: .bnet)"
    )
    parser.add_argument(
        "--csv", default=None, metavar="FILE",
        help="Write the summary table (model, count) to FILE as CSV"
    )
    args = parser.parse_args()

    try:
        import mpbn
    except ImportError:
        print("Error: mpbn is not installed. Run: pip install mpbn", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.folder):
        print(f"Error: '{args.folder}' is not a valid folder.", file=sys.stderr)
        sys.exit(1)

    files = sorted(f for f in os.listdir(args.folder) if f.endswith(args.ext))

    if not files:
        print(f"No '{args.ext}' files found in {args.folder}", file=sys.stderr)
        sys.exit(1)

    reachable_from = None
    if args.init_state:
        reachable_from = parse_state(args.init_state)

    results = []
    for fname in files:
        fpath = os.path.join(args.folder, fname)
        try:
            count = count_for_model(
                mpbn, fpath,
                fixed_points=args.fixed_points,
                reachable_from=reachable_from,
            )
            results.append({"file": fname, "n_attractors": count, "error": None})
            print(f"{fname:<40} {count}")
        except Exception as e:
            print(f"[WARNING] {fname}: {e}", file=sys.stderr)
            results.append({"file": fname, "n_attractors": None, "error": str(e)})

    df = pd.DataFrame(results)

    valid = df[df["n_attractors"].notna()]
    n_ok = len(valid)
    n_failed = len(df) - n_ok

    print("-" * 60)
    print(f"Models analyzed successfully : {n_ok}")
    if n_failed:
        print(f"Models failed                : {n_failed}")

    if n_ok > 0:
        counts = valid["n_attractors"]
        print(f"Attractor count -> min={counts.min()} max={counts.max()} "
              f"mean={counts.mean():.2f} median={counts.median()}")

        # Distribution of attractor counts across all models
        print("\nDistribution (n_attractors : number of models):")
        print(counts.value_counts().sort_index())

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nSummary exported to: {args.csv}")


if __name__ == "__main__":
    main()
