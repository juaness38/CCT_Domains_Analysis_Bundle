# cctdomain-peptide-networks

*[Versión en español](README.es.md)*

Residue interaction network analysis for CCT domain / peptide complexes sampled
by molecular dynamics.

The input to this pipeline is a molecular dynamics trajectory of a complex in
which one chain holds a CCT domain and another holds the peptide bound to it.
What the pipeline measures is which residue contacts hold that interface
together, and how reproducible those contacts are across independent samples of
the trajectory.

The scripts are protein-agnostic. They have been used on the CCT domains of
SPAK, WNK4, WNK2, NRBP1, OSR1 and FERRY3, and nothing about any of those
proteins is written into the code.

## How the pipeline works

Interaction networks computed from a single structure are not reproducible: a
contact may be present in one frame and absent in the next. The pipeline deals
with this by sampling the trajectory several times over and keeping only what
survives in the majority of samples.

1. **Sampling.** A set of frames spread evenly across the trajectory is split
   into several ensembles. The split is interleaved rather than block-wise, so
   each ensemble covers the whole trajectory instead of one stretch of it. This
   makes the ensembles independent samples of the same dynamics.
2. **Network calculation.** Each ensemble is turned into a residue interaction
   network. This step is done outside this repository, with a residue
   interaction network generator that writes Cytoscape JSON.
3. **Consensus.** The per-ensemble networks are merged. Every contact receives a
   consensus score: the fraction of networks in which it appears. A contact seen
   in all of them is reproducible; one seen in a single network is not.
4. **Filtering.** The consensus network is reduced to the contacts that cross
   the interface, above a frequency threshold.
5. **Renumbering.** Models number the CCT domain from 1, while the literature
   numbers it by its position in the full-length protein. This step restores the
   canonical numbering.
6. **Representative state.** From the ensemble, the single structure whose
   contact distances are closest to the consensus is selected, for figures and
   further analysis.

## Repository layout

```
scripts/                 the five analysis scripts
software_versions.tsv    inventory of every program used, with versions and evidence
```

The scripts are generalized implementations: the input paths, the chain roles,
the thresholds and the protein itself are command-line arguments rather than
constants written into the code.

They derive from earlier per-system implementations, one copy per protein, which
are not distributed here. Those earlier implementations produced the published
results, and the versions in this repository were checked against them on real
data and reproduce them exactly. See [Verification](#verification).

## The scripts

### `sample_ensembles.py`

Splits a trajectory into interleaved ensembles.

| | |
|---|---|
| Input | A trajectory readable by MDAnalysis, typically a multi-MODEL PDB. A separate topology file may be supplied when the format needs one. |
| Output | One multi-MODEL PDB per ensemble, plus a text file recording which frame indices went into each. |

Parameters control how many frames are selected in total, how many ensembles
they are divided into, the output prefix, and which atoms are written. Chain
identifiers present in the input are preserved, which matters because every
later step relies on them.

### `build_consensus.py`

Merges per-ensemble networks into one consensus network.

| | |
|---|---|
| Input | Two or more residue interaction networks in Cytoscape JSON form. |
| Output | The consensus network as JSON, plus an edge table and a node table as CSV. |

Edges are treated as undirected, so a contact recorded in either direction is
the same edge. Within one network, duplicate edges sharing a residue pair and an
interaction type are collapsed and their frequency and distance averaged. Across
networks, each edge receives a consensus score, along with the mean, minimum and
maximum of its frequency and the mean of its distance, and a record of which
input files it came from.

A minimum support parameter discards edges that appear in too few networks.

### `filter_consensus.py`

Reduces the consensus edge table to the interface contacts.

| | |
|---|---|
| Input | The edge table written by `build_consensus.py`. |
| Output | The same table with only the rows that pass the filters. |

Three conditions are applied together: the mean frequency must reach a
threshold, the interaction type must not be one of the excluded ones, and the
contact must run between the CCT domain chain and the peptide chain in either
direction. Contacts internal to one chain are dropped. The chains, the
threshold and the excluded interaction types are all parameters.

### `renumber_network.py`

Restores canonical numbering to the CCT domain.

| | |
|---|---|
| Input | The consensus edge table, the consensus network JSON, and the full-length FASTA of the CCT domain protein. |
| Output | A filtered and renumbered edge table, a renumbered network JSON, and a FASTA of the domain fragment. |

The CCT domain sequence is reconstructed from the network's own nodes, located
within the full-length FASTA, and every domain residue identifier is shifted by
the resulting offset. The peptide chain keeps the numbering it already has,
because it is a separate molecule whose local numbering is the meaningful one.
The same filters as `filter_consensus.py` are applied on the way through.

The script stops rather than guessing if the domain sequence is absent from the
FASTA, or if it occurs more than once and the mapping would be ambiguous.

The header of the fragment FASTA is derived from the header of the input FASTA,
so the fragment is always labelled with the protein it actually came from.

### `select_representative_state.py`

Selects the ensemble member closest to the consensus geometry.

| | |
|---|---|
| Input | A multi-MODEL PDB, a table of reference interactions as `.xlsx` or `.csv`, and the full-length FASTA. |
| Output | The selected structure as a single-model PDB with the CCT domain renumbered, and a report of every model's distances and scores. The whole ensemble can optionally be written out renumbered as well. |

Every listed interaction is measured in every model, and the model with the
smallest sum of squared deviations from the reference distances is selected.
Distance is defined by interaction type: between aromatic ring centroids for
PIPISTACK, between charged group centroids for IONIC, and as the shortest polar
heavy-atom pair for HBOND and anything else.

The reference table needs the columns `interaction`, `source`, `target` and
`average_distance_mean`.

## Conventions

### Chains

Throughout the pipeline, one chain holds the CCT domain and another holds the
peptide. Chains A and B respectively are the defaults, and both are parameters
wherever they are used.

Only the CCT domain is ever renumbered. The peptide keeps its local numbering.

### Residue identifiers

Residues are written as `CHAIN/NUMBER/RESNAME`, for instance a chain, a residue
number and a three-letter residue name separated by slashes. In the reference
table read by `select_representative_state.py`, the number may also be a range
written `START-END`, in which case the closest valid residue pair is used.

### Renumbering offsets

An offset of zero is a normal result, not a failure. It means the CCT domain
sits at the N-terminus of its protein, so the model numbering was already
canonical. Among the systems this pipeline has been used on, the offsets range
from zero to several hundred residues. The scripts report the offset they
derived and continue.

## Requirements

Python 3.8 or later, and:

| Package | Needed by |
|---|---|
| MDAnalysis | `sample_ensembles.py` |
| NumPy | `sample_ensembles.py`, `select_representative_state.py` |
| pandas | `select_representative_state.py` |
| openpyxl | `select_representative_state.py`, only to read `.xlsx` reference tables |

`build_consensus.py`, `filter_consensus.py` and `renumber_network.py` use only
the Python standard library.

The dependencies are listed in `requirements.txt`. Note that they are not all
needed at once: the two network steps run anywhere Python runs.

## Verification

The generalized scripts were checked against the published WNK4 CCT2 replica 2
analysis, by running them on the same inputs and comparing with the files the
original scripts produced.

| Script | Result |
|---|---|
| `build_consensus.py` | 605 edges and 100 nodes, matching the original, with no differing value in any numeric column |
| `filter_consensus.py` | the same 21 interface contacts as the original |
| `renumber_network.py` | the same 21 contacts, renumbered identically, with no differing cell; the renumbered network JSON identical in nodes and edges |
| `select_representative_state.py` | the same model selected, with the same sum of squared deviations and the same selected structure atom for atom |
| `sample_ensembles.py` | the frame selection and the interleaved split reproduce the recorded frame indices of all 22 archived runs, across every system. Writing the ensemble PDB files themselves is unverified, because MDAnalysis was unavailable in the environment where the check was run |

## License

MIT. See [LICENSE](LICENSE).
