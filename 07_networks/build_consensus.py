#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a consensus residue interaction network from several per-ensemble networks.

Each input is a residue interaction network in Cytoscape JSON form, with an
"elements" object holding "nodes" and "edges". Edges are treated as undirected,
so a contact between two residues is the same edge whichever way round it was
recorded.

Within a single network, duplicate edges sharing the same residue pair and the
same interaction type are collapsed and their frequency and distance averaged.
Across networks, every edge gets a consensus score: the fraction of input
networks in which it appears.

Writes the consensus network as JSON plus two CSV tables, one of edges and one
of nodes.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("networks", nargs="+",
                   help="Input network JSON files, one per ensemble")
    p.add_argument("-o", "--output-prefix", required=True,
                   help="Prefix for the three output files")
    p.add_argument("--min-support", type=int, default=1,
                   help="Keep only edges present in at least this many input "
                        "networks (default: 1, meaning keep everything)")
    return p.parse_args(argv)


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def canonical_edge_key(source, target, interaction):
    """Order the residue pair so the edge is undirected."""
    first, second = sorted([str(source), str(target)])
    return (first, second, str(interaction).upper())


def load_network(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    elements = data.get("elements", {})
    return elements.get("nodes", []), elements.get("edges", [])


NODE_FIELDS = ["id", "name", "chain", "resi", "resn",
               "degree", "x_coord", "y_coord", "z_coord", "dssp"]


def collect_node_attributes(known, nodes):
    """Remember each node's attributes the first time the node is seen."""
    for node in nodes:
        data = node.get("data", {})
        node_id = data.get("id") or data.get("name") or data.get("value")
        if not node_id:
            continue
        node_id = str(node_id)
        if node_id in known:
            continue
        entry = {field: data.get(field, "") for field in NODE_FIELDS}
        entry["id"] = node_id
        entry["name"] = data.get("name", node_id)
        known[node_id] = entry


def summarize_network_edges(edges):
    """Collapse duplicate edges inside one network, averaging their values."""
    gathered = defaultdict(lambda: {"frequencies": [], "distances": []})

    for edge in edges:
        data = edge.get("data", {})
        source = data.get("source")
        target = data.get("target")
        interaction = data.get("interaction")
        if not source or not target or not interaction:
            continue

        key = canonical_edge_key(source, target, interaction)
        frequency = safe_float(data.get("frequency"))
        distance = safe_float(data.get("distance"))
        if frequency is not None:
            gathered[key]["frequencies"].append(frequency)
        if distance is not None:
            gathered[key]["distances"].append(distance)

    summary = {}
    for key, values in gathered.items():
        frequencies = values["frequencies"]
        distances = values["distances"]
        summary[key] = {
            "frequency": sum(frequencies) / len(frequencies) if frequencies else None,
            "distance": sum(distances) / len(distances) if distances else None,
        }
    return summary


def build_consensus(paths, min_support=1):
    total_networks = len(paths)
    node_attributes = {}
    edge_observations = defaultdict(list)

    for path in paths:
        print(f"Reading: {path}")
        nodes, edges = load_network(path)
        collect_node_attributes(node_attributes, nodes)
        for key, stats in summarize_network_edges(edges).items():
            edge_observations[key].append({
                "file": os.path.basename(path),
                "frequency": stats["frequency"],
                "distance": stats["distance"],
            })

    consensus_edges = []
    for key, observations in edge_observations.items():
        source, target, interaction = key
        net_count = len(observations)
        if net_count < min_support:
            continue

        frequencies = [o["frequency"] for o in observations if o["frequency"] is not None]
        distances = [o["distance"] for o in observations if o["distance"] is not None]

        def rounded(value):
            return round(value, 6) if value is not None else ""

        label = f"{source} ({interaction}) {target}"
        consensus_edges.append({
            "consensus_score": round(net_count / total_networks, 6),
            "distance_mean": rounded(sum(distances) / len(distances) if distances else None),
            "frequency_max": rounded(max(frequencies) if frequencies else None),
            "frequency_mean": rounded(sum(frequencies) / len(frequencies) if frequencies else None),
            "frequency_min": rounded(min(frequencies) if frequencies else None),
            "interaction": interaction,
            "key": label,
            "name": label,
            "net_count": net_count,
            "present_in": ";".join(sorted(o["file"] for o in observations)),
            "selected": "TRUE",
            "shared interaction": interaction,
            "shared name": label,
            "source": source,
            "target": target,
        })

    consensus_edges.sort(
        key=lambda e: (-e["net_count"], e["interaction"], e["source"], e["target"])
    )

    used_nodes = {e[end] for e in consensus_edges for end in ("source", "target")}
    consensus_nodes = [
        {"data": node_attributes.get(node_id, {"id": node_id, "name": node_id})}
        for node_id in sorted(used_nodes)
    ]

    network = {
        "data": [],
        "directed": False,
        "multigraph": False,
        "elements": {
            "nodes": consensus_nodes,
            "edges": [
                {"data": dict(edge, id=f"e{i}")}
                for i, edge in enumerate(consensus_edges, start=1)
            ],
        },
    }
    return network, consensus_edges, consensus_nodes


EDGE_FIELDS = ["consensus_score", "distance_mean", "frequency_max", "frequency_mean",
               "frequency_min", "interaction", "key", "name", "net_count",
               "present_in", "selected", "shared interaction", "shared name",
               "source", "target"]


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    args = parse_args(argv)

    missing = [path for path in args.networks if not os.path.exists(path)]
    if missing:
        sys.exit("These input files were not found:\n" + "\n".join(missing))

    if args.min_support > len(args.networks):
        sys.exit(f"--min-support is {args.min_support} but only "
                 f"{len(args.networks)} networks were given; no edge could qualify.")

    network, edges, nodes = build_consensus(args.networks, args.min_support)

    json_out = f"{args.output_prefix}.json"
    edges_out = f"{args.output_prefix}_edges.csv"
    nodes_out = f"{args.output_prefix}_nodes.csv"

    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(network, handle, indent=2, ensure_ascii=False)
    write_csv(edges_out, EDGE_FIELDS, edges)
    write_csv(nodes_out, NODE_FIELDS, [node.get("data", {}) for node in nodes])

    print("\nDone.")
    print(f"Networks read: {len(args.networks)}")
    print(f"JSON:  {json_out}")
    print(f"EDGES: {edges_out}")
    print(f"NODES: {nodes_out}")
    print(f"Nodes: {len(nodes)}")
    print(f"Edges: {len(edges)}")


if __name__ == "__main__":
    main()
