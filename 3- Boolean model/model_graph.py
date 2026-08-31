#!/usr/bin/env python3

import argparse
import os
import bonesistools as bt

parser = argparse.ArgumentParser(
    description="Generate ensemble influence graph from BoNesis models"
)

parser.add_argument(
    "--results",
    required=True,
    help="Results folder (e.g. results_cite-seq_i2)"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output PDF filename"
)

args = parser.parse_args()

model_dir = f"{args.results}/infer/bn/min/model"
model_out = f"{args.results}/{args.output}"


print(f"Loading models from: {model_dir}")

bns = bt.bpy.bn.read_bnet_directory(
    model_dir,
    recursive=True
)

dot = bns.to_pydot(remove_isolated_nodes=True)

success = dot.write(model_out, format="pdf")

if success:
    print(f"PDF written to: {args.output}")
else:
    print("Failed to generate PDF")
