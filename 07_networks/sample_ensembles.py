#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Split a molecular dynamics trajectory into interleaved structural ensembles.

Reads a multi-frame trajectory, picks a number of uniformly spaced frames and
distributes them round-robin into several multi-MODEL PDB files. Interleaving
means consecutive selected frames land in different ensembles, so the ensembles
sample the same stretch of the trajectory instead of consecutive blocks.

Chain identifiers present in the input are preserved in the output.
"""

import argparse
import sys

import numpy as np
import MDAnalysis as mda


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("trajectory",
                   help="Input trajectory readable by MDAnalysis (multi-MODEL PDB, "
                        "or a topology/trajectory pair via --topology)")
    p.add_argument("--topology", default=None,
                   help="Separate topology file, when the trajectory format needs one")
    p.add_argument("-o", "--output-prefix", default=None,
                   help="Prefix for the generated files (default: the trajectory "
                        "file name without its extension)")
    p.add_argument("-n", "--frames", type=int, default=300,
                   help="Total number of frames to select (default: 300)")
    p.add_argument("-e", "--ensembles", type=int, default=3,
                   help="Number of ensembles to split them into (default: 3)")
    p.add_argument("--selection", default="all",
                   help="MDAnalysis atom selection to write (default: all)")
    return p.parse_args(argv)


def uniform_positions(n_frames, n_select):
    """Indices of n_select frames spread evenly over n_frames, without repeats."""
    positions = np.unique(np.rint(np.linspace(0, n_frames - 1, n_select)).astype(int))
    if len(positions) < n_select:
        # Rounding collapsed neighbouring picks; fall back to flooring over the
        # half-open interval, which spreads them slightly differently.
        positions = np.unique(
            np.floor(np.linspace(0, n_frames, n_select, endpoint=False)).astype(int)
        )
    return positions


def main(argv=None):
    args = parse_args(argv)

    if args.frames < 1 or args.ensembles < 1:
        sys.exit("--frames and --ensembles must both be at least 1.")

    universe = (mda.Universe(args.topology, args.trajectory)
                if args.topology else mda.Universe(args.trajectory))
    n_frames = len(universe.trajectory)

    prefix = args.output_prefix
    if prefix is None:
        prefix = args.trajectory.rsplit(".", 1)[0]

    print(f"Loaded: {args.trajectory}")
    print(f"Frames found: {n_frames}")

    if n_frames < args.frames:
        sys.exit(f"The trajectory holds {n_frames} frames but {args.frames} were "
                 f"requested. Lower --frames, or use a longer trajectory.")

    positions = uniform_positions(n_frames, args.frames)
    if len(positions) < args.frames:
        sys.exit(f"Could not select {args.frames} distinct frames out of {n_frames}; "
                 f"only {len(positions)} are distinguishable.")

    subsets = [positions[i::args.ensembles] for i in range(args.ensembles)]
    sizes = sorted({len(s) for s in subsets})
    if len(sizes) > 1:
        print(f"Note: {args.frames} frames do not divide evenly into "
              f"{args.ensembles} ensembles; sizes range from {sizes[0]} to {sizes[-1]}.")

    atoms = universe.select_atoms(args.selection)
    if atoms.n_atoms == 0:
        sys.exit(f"The selection {args.selection!r} matched no atoms.")

    for i, subset in enumerate(subsets, start=1):
        outname = f"{prefix}_ensemble_{i}.pdb"
        with mda.Writer(outname, multiframe=True, n_atoms=atoms.n_atoms) as writer:
            for frame_index in subset:
                universe.trajectory[frame_index]
                writer.write(atoms)
        print(f"Wrote: {outname} ({len(subset)} frames)")

    index_file = f"{prefix}_selected_frames.txt"
    with open(index_file, "w", encoding="utf-8") as handle:
        handle.write(f"# {len(positions)} uniformly spaced frame indices\n")
        handle.write(",".join(map(str, positions)) + "\n\n")
        for i, subset in enumerate(subsets, start=1):
            handle.write(f"# ensemble_{i}\n")
            handle.write(",".join(map(str, subset)) + "\n\n")

    print(f"Frame indices written to: {index_file}")


if __name__ == "__main__":
    main()
