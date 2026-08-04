#!/usr/bin/env python3

# USER 
# 1) conda activate bonesis 
# 2) For utilisation: python plot_models.py --input models_selected --output models_selected/graph


import argparse
from pathlib import Path

import bonesistools as bt

GRAPH_ATTR = {
    "ratio": "compress",
    "overlap": "prism",
    "sep": "+0",
    "esep": "+0",
    "K": "0.35",
    "ranksep": "0.6",
    "pack": "true",
    "rankdir": "TB",
    "splines": "curve",
}

BASE_OPTIONS = {
    "edge_label": "frequency",
    "node_style": "stability",
    "edge_style": "frequency",
    "preserve_feedback": True,
    "include_selfloops": False,
    "min_frequency": 0,
    "program": "dot",
    "graph_attr": GRAPH_ATTR,
    "node_attr": {"fontsize": "20"},
    "edge_attr": {"fontsize": "20"},
}


def export(graph, outfile, collapse=None, drop_isolates=True):

    dot = graph.to_pydot(
        **BASE_OPTIONS,
        collapse=collapse,
        drop_isolates=drop_isolates,
    )

    dot.write_jpg(str(outfile))
    print(f"Saved {outfile}")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="Directory containing .bnet files")
    parser.add_argument("-o", "--output", default="influence_graph")
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(exist_ok=True)

    bns = bt.logic.io.read_bnet_directory(args.input, recursive=True)
    graph = bns.to_influence_graph()

    export(
        graph,
        outdir / "aggregate.jpg",
        collapse=None,
        drop_isolates=True,
    )

    export(
        graph,
        outdir / "aggregate_with_isolates.jpg",
        collapse=None,
        drop_isolates=False,
    )

    export(
        graph,
        outdir / "feedback_core.jpg",
        collapse="feedback",
        drop_isolates=True,
    )

    export(
        graph,
        outdir / "function_families.jpg",
        collapse="family",
        drop_isolates=True,
    )


if __name__ == "__main__":
    main()
