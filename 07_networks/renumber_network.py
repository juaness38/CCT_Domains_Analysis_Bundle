#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Renumber the CCT domain of a consensus network to canonical sequence numbering.

Models of a CCT domain usually number their residues from 1, while the literature
refers to the residues by their position in the full-length protein. This script
restores the canonical numbering.

It reconstructs the CCT domain sequence from the network's own nodes, locates it
inside the full-length FASTA and shifts every CCT domain residue identifier by the
resulting offset. The peptide chain keeps the numbering it already has, because it
is a separate molecule whose local numbering is the meaningful one.

An offset of zero is a normal outcome, not an error: it means the domain sits at
the N-terminus of the protein, so the model numbering was canonical already.

The same filters as filter_consensus.py are applied to the edge table on the way
through. Three files are written: the filtered and renumbered edge table, the
renumbered network JSON, and a FASTA of the domain fragment, whose header is
derived from the input FASTA so it always names the right protein.
"""

import argparse
import csv
import json
import os
import re
import sys
from copy import deepcopy

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

RENUMBERED_TEXT_FIELDS = ("key", "name", "shared name")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("edges_csv", help="Edge table from build_consensus.py")
    p.add_argument("network_json", help="Network JSON from build_consensus.py")
    p.add_argument("fasta", help="Full-length FASTA of the CCT domain protein")
    p.add_argument("--cct-domain-chain", default="A",
                   help="Chain holding the CCT domain, the one renumbered (default: A)")
    p.add_argument("--peptide-chain", default="B",
                   help="Chain holding the peptide, left untouched (default: B)")
    p.add_argument("--min-frequency", type=float, default=0.1,
                   help="Minimum frequency_mean to keep a row (default: 0.1)")
    p.add_argument("--exclude-interaction", nargs="*", default=["VDW"], metavar="TYPE",
                   help="Interaction types to drop (default: VDW)")
    p.add_argument("-o", "--output-prefix", default=None,
                   help="Prefix for the output files (default: derived from the inputs)")
    return p.parse_args(argv)


def read_fasta(path):
    """Return (header, sequence) of the first record, sequence upper-cased."""
    header = ""
    parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:]
            else:
                parts.append(line.replace(" ", "").upper())
    sequence = "".join(parts)
    if not sequence:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return header, sequence


def fasta_label(header):
    """A short identifier for the protein, taken from its FASTA header.

    Handles the UniProt form sp|ACCESSION|NAME, and otherwise falls back to the
    first whitespace-separated token. This is what keeps the fragment FASTA from
    being labelled with whatever protein the script was first written for.
    """
    header = (header or "").strip()
    if not header:
        return "protein"
    if header.startswith(("sp|", "tr|")) and header.count("|") >= 2:
        return header.split("|")[2].split()[0]
    return header.split()[0]


def parse_residue_id(value):
    """CHAIN/NUMBER/RESNAME -> (chain, number, resname); number is None if absent."""
    parts = str(value).strip().split("/")
    if len(parts) < 3:
        return "", None, ""
    try:
        number = int(parts[1])
    except (TypeError, ValueError):
        number = None
    return parts[0].upper(), number, parts[2].upper()


def chain_sequence_from_network(data, chain):
    """Rebuild a chain's sequence from the network nodes.

    Returns (sequence, first_residue_number, last_residue_number).
    """
    residues = {}
    for node in data.get("elements", {}).get("nodes", []):
        node_data = node.get("data", {})
        if str(node_data.get("chain", "")).upper() != chain:
            continue
        try:
            number = int(node_data.get("resi"))
        except (TypeError, ValueError):
            continue
        resname = str(node_data.get("resn", "")).upper()
        if resname not in AA3_TO_AA1:
            raise ValueError(f"Unrecognised residue in chain {chain}: "
                             f"{resname} at position {number}")
        residues[number] = resname

    if not residues:
        raise ValueError(f"No chain {chain} residues found in the network JSON. "
                         f"Check --cct-domain-chain.")

    positions = sorted(residues)
    if positions != list(range(positions[0], positions[-1] + 1)):
        raise ValueError(f"Chain {chain} is not numbered contiguously in the network "
                         f"JSON, so the canonical numbering cannot be inferred "
                         f"unambiguously.")

    sequence = "".join(AA3_TO_AA1[residues[p]] for p in positions)
    return sequence, positions[0], positions[-1]


def locate_in_fasta(fragment, full_sequence):
    """Find the fragment in the full sequence; return 1-based start and end."""
    first = full_sequence.find(fragment)
    if first < 0:
        raise ValueError("The CCT domain sequence from the network was not found "
                         "in the FASTA. The script stops rather than guessing the "
                         "numbering; check that the FASTA is the right protein.")
    if full_sequence.find(fragment, first + 1) >= 0:
        raise ValueError("The CCT domain sequence occurs more than once in the "
                         "FASTA, so the canonical numbering is ambiguous.")
    return first + 1, first + len(fragment)


def make_renumberer(chain, offset):
    """Return functions that shift chain residue numbers by a fixed offset."""
    pattern = re.compile(rf"\b{re.escape(chain)}/(\d+)/")

    def renumber_id(value):
        value_chain, number, resname = parse_residue_id(value)
        if value_chain != chain or number is None:
            return value
        return f"{chain}/{number + offset}/{resname}"

    def renumber_text(text):
        if not isinstance(text, str):
            return text
        return pattern.sub(lambda m: f"{chain}/{int(m.group(1)) + offset}/", text)

    return renumber_id, renumber_text


def renumber_network(data, chain, offset, renumber_id, renumber_text):
    out = deepcopy(data)
    elements = out.get("elements", {})

    for node in elements.get("nodes", []):
        node_data = node.get("data", {})
        if str(node_data.get("chain", "")).upper() != chain:
            continue
        try:
            number = int(node_data.get("resi"))
        except (TypeError, ValueError):
            continue
        new_number = number + offset
        resname = str(node_data.get("resn", "")).upper()
        node_data["resi"] = new_number
        for field in ("id", "name"):
            if field in node_data:
                node_data[field] = f"{chain}/{new_number}/{resname}"

    for edge in elements.get("edges", []):
        edge_data = edge.get("data", {})
        for field in ("source", "target"):
            if field in edge_data:
                edge_data[field] = renumber_id(edge_data[field])
        for field in RENUMBERED_TEXT_FIELDS:
            if field in edge_data:
                edge_data[field] = renumber_text(edge_data[field])

    return out


def renumber_edge_table(input_csv, output_csv, pair, min_frequency, excluded,
                        renumber_id, renumber_text):
    with open(input_csv, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing = {"frequency_mean", "interaction", "source", "target"} - set(fieldnames)
        if missing:
            raise ValueError("Missing required columns in the edge table: " +
                             ", ".join(sorted(missing)))
        rows = list(reader)

    kept = []
    for row in rows:
        try:
            frequency = float(row.get("frequency_mean", ""))
        except (TypeError, ValueError):
            continue
        if frequency < min_frequency:
            continue
        if str(row.get("interaction", "")).strip().upper() in excluded:
            continue
        chains = {parse_residue_id(row.get("source"))[0],
                  parse_residue_id(row.get("target"))[0]}
        if chains != pair:
            continue

        new_row = dict(row)
        for field in ("source", "target"):
            new_row[field] = renumber_id(row.get(field, ""))
        for field in RENUMBERED_TEXT_FIELDS:
            if field in new_row:
                new_row[field] = renumber_text(new_row[field])
        kept.append(new_row)

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(kept)

    return len(rows), len(kept)


def write_fragment_fasta(path, label, chain, sequence, start, end, width=60):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{label}_CCT_domain_chain_{chain}_{start}-{end}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i:i + width] + "\n")


def main(argv=None):
    args = parse_args(argv)

    for path in (args.edges_csv, args.network_json, args.fasta):
        if not os.path.exists(path):
            sys.exit(f"File not found: {path}")

    domain = args.cct_domain_chain.strip().upper()
    peptide = args.peptide_chain.strip().upper()
    if domain == peptide:
        sys.exit(f"--cct-domain-chain and --peptide-chain must differ; both are {domain!r}.")
    excluded = {t.strip().upper() for t in args.exclude_interaction}

    header, full_sequence = read_fasta(args.fasta)
    with open(args.network_json, "r", encoding="utf-8") as handle:
        network = json.load(handle)

    try:
        sequence, local_start, local_end = chain_sequence_from_network(network, domain)
        canonical_start, canonical_end = locate_in_fasta(sequence, full_sequence)
    except ValueError as error:
        sys.exit(str(error))

    offset = canonical_start - local_start
    label = fasta_label(header)

    print("Mapping")
    print(f"  FASTA:            {header}")
    print(f"  protein label:    {label}")
    print(f"  chain {domain} local:     {local_start}-{local_end} ({len(sequence)} aa)")
    print(f"  canonical:        {canonical_start}-{canonical_end}")
    print(f"  offset:           {offset:+d}")
    if offset == 0:
        print("  The domain already carries canonical numbering, so renumbering "
              "changes no residue number. This is expected when the domain sits "
              "at the N-terminus of the protein.")

    prefix = args.output_prefix
    if prefix is None:
        prefix = os.path.splitext(args.edges_csv)[0]
    edges_out = f"{prefix}_filtered_renumbered.csv"
    json_out = f"{os.path.splitext(args.network_json)[0]}_renumbered.json"
    fasta_out = f"{prefix}_{label}_CCT_domain_{canonical_start}_{canonical_end}.fasta"

    renumber_id, renumber_text = make_renumberer(domain, offset)

    try:
        total_rows, kept_rows = renumber_edge_table(
            args.edges_csv, edges_out, {domain, peptide},
            args.min_frequency, excluded, renumber_id, renumber_text)
    except ValueError as error:
        sys.exit(str(error))

    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(renumber_network(network, domain, offset, renumber_id, renumber_text),
                  handle, indent=2, ensure_ascii=False)

    write_fragment_fasta(fasta_out, label, domain, sequence, canonical_start, canonical_end)

    print("\nDone. Files written:")
    print(f"  filtered and renumbered edges: {edges_out}")
    print(f"  renumbered network:            {json_out}")
    print(f"  CCT domain fragment FASTA:     {fasta_out}")
    print(f"\nEdge table: {total_rows} rows in, {kept_rows} kept")
    print(f"Filters: frequency_mean >= {args.min_frequency}; "
          f"excluded {', '.join(sorted(excluded)) if excluded else 'nothing'}; "
          f"chains {domain}-{peptide} only")
    if not kept_rows:
        print("Warning: no row survived the filters. Check the chain identifiers.")


if __name__ == "__main__":
    main()
