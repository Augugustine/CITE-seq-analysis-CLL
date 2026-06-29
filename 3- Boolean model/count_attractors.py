#!/usr/bin/env python3
"""
count_attractors.py - Count the attractors of a Boolean network (BoolNet format)

Usage:
    python count_attractors.py <file.bnet>
    python count_attractors.py <file.bnet> --fixed-points
    python count_attractors.py <file.bnet> --from "a=0,b=1,c=0"
    python count_attractors.py <file.bnet> --output attractors.txt

Options:
    --fixed-points      Count fixed points only
    --output FILE       Write attractor vectors to FILE
    --from STATE        Restrict to attractors reachable from STATE
                        Format: "node=value,node=value,..."
"""

import sys
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

def main():
    parser = argparse.ArgumentParser(
        description="Count attractors of a Boolean network (.bnet) using mpbn"
    )
    parser.add_argument("bnet_file", help="Boolean network file in BoolNet format (.bnet)")
    parser.add_argument(
        "--fixed-points", action="store_true",
        help="Count fixed points only"
    )
    parser.add_argument(
        "--output", "-o", default=None, metavar="FILE",
        help="Write attractor vectors to FILE (e.g. attractors.txt)"
    )
    parser.add_argument(
        "--from", dest="init_state", default=None, metavar="STATE",
        help="Initial state, e.g.: 'a=0,b=1,c=0'"
    )
    args = parser.parse_args()

    try:
        import mpbn
    except ImportError:
        print("Error: mpbn is not installed. Run: pip install mpbn", file=sys.stderr)
        sys.exit(1)

    try:
        mbn = mpbn.MPBooleanNetwork(args.bnet_file)
    except FileNotFoundError:
        print(f"Error: file not found: {args.bnet_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading network: {e}", file=sys.stderr)
        sys.exit(1)

    reachable_from = None
    if args.init_state:
        reachable_from = parse_state(args.init_state)

    kwargs = {}
    if reachable_from:
        kwargs["reachable_from"] = reachable_from

    if args.fixed_points:
        count = mbn.count_fixedpoints(**kwargs)
        iterator = mbn.fixedpoints(**kwargs)
    else:
        count = mbn.count_attractors(**kwargs)
        iterator = mbn.attractors(**kwargs)

    # Main output: just the number
    print(count)

    # Write vectors to file if requested
    # Write attractors as CSV
    if args.output:
    	try:
        	# Dictionary
        	attractors = list(iterator)

        	# DataFrame
        	df = pd.DataFrame(attractors)

        	df = df[sorted(df.columns)]

        	# Row names
        	df.index = [f"Attractor{i+1}" for i in range(len(df))]
        	df.index.name = "Attractor"

        	# Save
        	df.to_csv(args.output)

        	print(f"CSV written to: {args.output}", file=sys.stderr)

    	except Exception as e:
        	print(f"Error writing output file: {e}", file=sys.stderr)
        	sys.exit(1)


if __name__ == "__main__":
    main()
    
