"""Compute UniProt residue offsets for 6 MD systems via pairwise local alignment.

DEPRECATED / HISTORICAL: this script produced uniprot_offset_map.json, whose
offsets are superseded by 03_data/uniprot_mapping.json. The UniProt sequences
embedded below are NOT canonical for O95747 (1-aa insertion), Q96J92 (1-aa
deletion) and Q9Y3S1 (different 841-aa sequence), which made the osr1, wnk4
and wnk2 offsets wrong. Do not run; kept only as provenance of the original
numbering. See _apply_uniprot_numbering.py for the canonical mapping.
"""
from __future__ import annotations

import json
import pathlib
from typing import Dict, List, Tuple

from Bio import Align

# ── constants ────────────────────────────────────────────────────────
BUNDLE_ROOT = pathlib.Path(
    r"C:\Users\busta\Downloads\MICA\_replot_md_controls_20260327\organized_by_system_bundle_20260327"
)
OUT_JSON = pathlib.Path(
    r"C:\Users\busta\Downloads\MICA\_replot_md_controls_20260327\uniprot_offset_map.json"
)

AA3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}

# ── UniProt sequences (verbatim, spaces/newlines stripped at runtime) ─
UNIPROT_DATA: Dict[str, Tuple[str, str, str]] = {
    # slug: (accession, gene, raw_sequence)
    "ferry_cctl": (
        "Q9NQ89", "FERRY3",
        "MKKNRERFCNREREFVYKFKVGSQCLELRVPLKFPVQENASHLHGRLMLLHSLPCFIEKDLKEALTQFIEEESLSDYDRD"
        "AEASLAAVKSGEVDLHQLASTWAKAYAETTLEHARPEEPSWDEDFADVYHDLIHSPASETLLNLEHNYFVSISELIGERD"
        "VELKKLRERQGIEMEKVMQELGKSLTDQDVNSLAAQHFESQQDLENKWSNELKQSTAIQKQEYQEWVIKLHQDLKNPNN"
        "SSLSEEIKVQPSQFRESVEAIGRIYEEQRKLEESFTIHLGAQLKTMHNLRLLRADMLDFCKHKRNHRSGVKLHRLQTALS"
        "LYSTSLCGLVLLVDNRINSYSGIKRDFATVCQECTDFHFPRIEEQLEVVQQVVLYARTQRRSKLKESLDSGNQNGGNDD"
        "KTKNAERNYLNVLPGEFYITRHSNLSEIHVAFHLCVDDHVKSGNITARDPAIMGLRNILKVCCTHDITTISIPLLLVHD"
        "MSEEMTIPWCLRRAELVFKCVKGFMMEMASWDGGISRTVQFLVPQSISEEMFYQLSNMLPQIFRVSSTLTLTSKH",
    ),
    "spak_cct": (
        "Q9UEW8", "STK39",
        "MAEPSGSPVHVQLPQQAAPVTAAAAAAPAAATAAPAPAAPAAPAPAPAPAAQAVGWPICRDAYELQEVIGSGATAVVQAA"
        "LCKPRQERVAIKRINLEKCQTSMDELLKEIQAMSQCSHPNVVTYYTSFVVKDELWLVMKLLSGGSMLDIIKYIVNRGEH"
        "KNGVLEEAIIATILKEVLEGLDYLHRNGQIHRDLKAGNILLGEDGSVQIADFGVSAFLATGGDVTRNKVRKTFVGTPCW"
        "MAPEVMEQVRGYDFKADMWSFGITAIELATGAAPYHKYPPMKVLMLTLQNDPPTLETGVEDKEMMKKYGKSFRKLLSLC"
        "LQKDPSKRPTAAELLKCKFFQKAKNREYLIEKLLTRTPDIAQRAKKVRRVPGSSGHLHKTEDGDWEWSDDEMDEKSEE"
        "GKAAFSQEKSRRVKEENPEIAVSASTIPEQIQSLSVHDSQGPPNANEDYREASSCAVNLVLRLRNSRKELNDIRFEFTPG"
        "RDTADGVSQELFSAGLVDGHDVVIVAANLQKIVDDPKALKTLTFKLASGCDGSEIPDEVKLIGFAQLSVS",
    ),
    "osr1_cct": (
        "O95747", "OXSR1",
        "MSEDSSALPWSINRDDYELQEVIGSGATAVVQAAYCAPKKEKVAIKRINLEKCQTSMDELLKEIQAMSQCHHPNIVSYYT"
        "SFVVKDELWLVMKLLSGGSVLDIIKHIVAKGEHKSGVLDESSTIATILREVLEGLEYLHKNGQIHRDVKAGNILLGEDGS"
        "VQIADFGVSAFLATGGDITRNKVRKTFVGTPCWMAPEVMEQVRGYDFKADIWSFGITAIELATGAAPYHKYPPMKVLML"
        "TLQNDPPSLETGVQDKEMLKKYGKSFRKMISLCLQKDPEKRPTAAELLRHKFFQKAKNKEFLQEKTLQRAPTISERAKKV"
        "RRVPGSSGRLHKTEDGGWEWSDDEFDEESEEGKAAISQLRSPRVKESISNSELFPTTDPVGTLLQVPEQISAHLPQPAGQ"
        "IATQPTQVSLPPTAEPAKTAQALSSGSGSQETKIPISLVLRLRNSKKELNDIRFEFTPGRDTAEGVSQELISAGLVDGRD"
        "LVIVAANLQKIVEEPQSNRSVTFKLASGVEGSDIPDDGKLIGFAQLSIS",
    ),
    "wnk4_cctl2": (
        "Q96J92", "WNK4",
        "MLASPATETTVLMSQTEADLALRPPPPLGTAGQPRLGPPPRRARRFSGKAEPRPRSSRLSRRSSVDLGLLSSWSLPASAP"
        "DPPDPPDSAGPGPARSPPPSSKEPPEGTWTEGAPVKAAEDSARPELPDSAVGPGSREPLRVPEAVALERRREQEEKEDMET"
        "QAVATSPDGRYLKFDIEIGNGSFKTVYRGLDTDTTVEVAWCELQTRKLSRAERQRFSEEVEMLKGLQHPNIVRFYDFWKS"
        "VLRGQVCIVLVTELMTSGTLKTYLRRFREMKPRVLQRWSRQILRGLHFLHSRVPPILHRDLKCDNVFITGPTGSVKIGD"
        "LGLATLKRASFAKSVIGTPEFMAPEMYEEKYDEAVDVYAFGMCMLEMATSEYPYSECQNAAQIYRKVTSGRKPNSFHKVK"
        "IPEVKEIIEGCIRTDKNERFTIQDLLAHAFFREERGVHVELAEEDDGEKPGLKLWLRMEDARRGGRPRDNQAIEFLFQLGR"
        "DAAEEVAQEMVALGLVCEADYQPVARAVRERVAAIQRKREKLRKARELEALPPEPGPPPATVPMAPGPPSVFPPEPEEPE"
        "ADQHQPFLFRHASYSSTTSDCETDGYLSSSGFLDASDPALQPPGGVPSSLAESHLCLPSAFALSIPRSGPGSDFSPGDSY"
        "ASDAASGLSDVGEGMGQMRRPPGRNLRRRPRSRLRVTSVSDQNDRVVECQLQTHNSKMVTFRFDLDGDSPEEIAAAMVY"
        "NEFILPSERDGFLRRIREIIQRVETLLKRDTGPMEAAEDTLSPQEEPAPLPALPVPLPDPSNEELQSSTSLEHRSWTAFS"
        "TSSSSPGTPLSPGNPFSPGTPISPGPIFPITSPPCHPSPSPFSPISSQVSSNPSPHPTSSPLPFSSSTPEFPVPLSQCPW"
        "SSLPTTSPPTFSPTCSQVTLSSPFFPPCPSTSSFPSTTAAPLLSLASAFSLAVMTVAQSLLSPSPGLLSQSPPAPPSPLP"
        "SLPLPPPVAPGGQESPSPHTAEVESEASPPPARPLPGEARLAPISEEGKPQLVGRFQVTSSKEPAEPLPLQPTSPTLSGSP"
        "KPSTPQLTSESSDTEDSAGGGPETREALAESDRAAEGLGAGVEEEGDDGKEPQVGGSPQPLSHPSPVWMNYSYSSLCLSS"
        "EESESSGEDEEFWAELQSLRQKHLSEVETLQTLQKKEIEDLYSRLGKQPPPGIVAPAAMLSSRQRRLSKGSFPTSRRNSL"
        "QRSEPPGPGIMRRNSLSGSSTGSQEQRASKGVTFAGDVGRM",
    ),
    "wnk2_cctl1": (
        "Q9Y3S1", "WNK2",
        "MDGDGGRRDVPGTLMEPGRGAGPAGMAEPRAKAARPGPQRFLRRSVVESDQEEPPGLEAAEAPGPQPPQPLQRRVLLLCK"
        "TRRLIAERARGRPAAPAPAALVAQPGAPGAPADAGPEPVGTQEPGPDPIAAAVETAPAPDGGPREEAAATVRKEDEGAAEA"
        "KPEPGRTRRDEPEEEEEDDEDDLKAVATSLDGRFLKFDIELGRGSFKTVYKGLDTETWVEVAWCELQDRKLTKLERQRFKE"
        "EAEMLKGLQHPNIVRFYDFWESSAKGKRCIVLVTELMTSGTLKTYLKRFKVMKPKVLRSWCRQILKGLLFLHTRTPPIIH"
        "RDLKCDNIFITGPTGSVKIGDLGLATLKRASFAKSVIGTPEFMAPEMYEEHYDESVDVYAFGMCMLEMATSEYPYSECQNA"
        "AQIYRKVTCGIKPASFEKVHDPEIKEIIGECICKNKEERYEIKDLLSHAFFAEDTGVRVELAEEDHGRKSTIALRLWVEDP"
        "KKLKGKPKDNGAIEFTFDLEKETPDEVAQEMIESGFFHESDVKIVAKSIRDRVALIQWRRERIWPALQPKEQQDVGSPDKA"
        "RGPPVPLQVQVTYHAQAGQPGPPEPEEPEADQHLLPPTLPTSATSLASDSTFDSGQGSTVYSDSQSSQQSVMLGSLADAAP"
        "SPAQCVCSPPVSEGPVLPQSLPSLGAYQQPTAAPGLPVGSVPAPACPPSLQQHFPDPAMSFAPVLPPPSTPMPTGPGQPA"
        "PPGQQPPPLAQPTPLPQVLAPQPVVPLQPVPPHLPPYLAPASQVGAPAQLKPLQMPQAPLQPLAQVPPQMPPIPVVPPIT"
        "PLAGIDGLPPALPDLPTATVPPVPPPQYFSPAVIL",
    ),
    "nrbp1_cct": (
        "Q9UHY1", "NRBP1",
        "MSEGESQTVLSSGSDPKVESSSSAPGLTSVSPPVTSTTSAASPEEEEESEDESEILEESP"
        "CGRWQKRREEVNQRNVPGIDSAYLAMDTEEGVEVVWNEVQFSERKNYKLQEEKVRAVFDNLIQLEHLNIVKFHKYWADIKE"
        "NKARVIFITEYMSSGSLKQFLKKTKKNHKTMNEKAWKRWCTQILSALSYLHSCDPPIIHGNLTCDTIFIQHNGLIKIGSVA"
        "PDTINNHVKTCREEQKNLHFFAPEYGEVTNVTTAVDIYSFGMCALEMAVLEIQGNGESSYVPQEAISSAIQLLEDPLQREF"
        "IQKCLQSEPARRPTARELLFHPALFEVPSLKLLAAHCIVGHQHMIPENALEEITKNMDTSAVLAEIPAGPGREPVQTLYSQS"
        "PALELDKFLEDVRNGIYPLTAFGLPRPQQPQQEEVTSPVVPPSVKTPTPEPAEVETRKVVLMQCNIESVEEGVKHHLTLLL"
        "KLEDKLNRHLSCDLMPNENIPELAAELVQLGFISEADQSRLTSLLEETLNKFNFARNSTLNSAAVTVSS",
    ),
}


def extract_ca_residues(pdb_path: pathlib.Path) -> List[Tuple[int, str]]:
    """Return list of (resSeq, resName_3letter) from CA ATOM records, chain A only."""
    seen: set = set()
    residues: List[Tuple[int, str]] = []
    with open(str(pdb_path), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21]
            if chain not in ("A", " "):
                continue
            res_name = line[17:20].strip()
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            key = (res_seq, chain)
            if key in seen:
                continue
            seen.add(key)
            residues.append((res_seq, res_name))
    residues.sort(key=lambda x: x[0])
    return residues


def build_sequence(residues: List[Tuple[int, str]]) -> str:
    return "".join(AA3TO1.get(res_name.upper(), "X") for _, res_name in residues)


def align_pdb_to_uniprot(
    pdb_seq: str,
    uniprot_seq: str,
) -> Tuple[float, int, int, int, float]:
    """Local pairwise alignment (PDB fragment -> full UniProt).

    Returns (score, query_start_0based, target_start_0based, target_end_0based, pct_coverage).
    target = UniProt, query = PDB.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(uniprot_seq, pdb_seq)
    best = next(iter(alignments))
    target_blocks = best.aligned[0]   # UniProt blocks
    query_blocks = best.aligned[1]    # PDB blocks
    target_start = int(target_blocks[0][0])
    target_end = int(target_blocks[-1][1])
    query_start = int(query_blocks[0][0])
    aligned_pdb_len = sum(b[1] - b[0] for b in query_blocks)
    pct_coverage = aligned_pdb_len / len(pdb_seq) * 100
    return float(best.score), query_start, target_start, target_end, pct_coverage


def main() -> None:
    results: Dict[str, dict] = {}

    for slug, (accession, gene, raw_seq) in UNIPROT_DATA.items():
        uniprot_seq = raw_seq.replace(" ", "").replace("\n", "")
        pdb_path = BUNDLE_ROOT / slug / "R1" / "dssp" / "last_frame_chainA.pdb"
        if not pdb_path.exists():
            pdb_path = BUNDLE_ROOT / slug / "R1" / "dssp" / "last_frame.pdb"

        print(f"\n{'=' * 64}")
        print(f"System : {slug}  ({accession} / {gene})")
        print(f"PDB    : {pdb_path}")

        ca_residues = extract_ca_residues(pdb_path)
        pdb_seq = build_sequence(ca_residues)
        first_res, last_res = ca_residues[0][0], ca_residues[-1][0]
        print(f"PDB CA : {len(ca_residues)} residues, resSeq {first_res}–{last_res}")
        print(f"PDB seq: {pdb_seq[:30]}...")
        print(f"UniProt: {uniprot_seq[:30]}...  ({len(uniprot_seq)} aa)")

        score, q_start, t_start, t_end, coverage = align_pdb_to_uniprot(pdb_seq, uniprot_seq)

        pdb_resSeq_start = ca_residues[q_start][0]
        uniprot_start_1based = t_start + 1
        uniprot_end_1based = t_end  # t_end is 0-based exclusive == 1-based inclusive

        offset = uniprot_start_1based - pdb_resSeq_start

        print(f"Score  : {score:.1f}   Coverage: {coverage:.1f}%")
        print(f"PDB seq idx {q_start} → resSeq {pdb_resSeq_start}")
        print(f"UniProt range: {uniprot_start_1based}–{uniprot_end_1based}")
        print(f"Offset (uniprot_pos - pdb_resSeq): {offset}")
        print("First 10 aligned positions (PDB resSeq → UniProt pos):")
        for i in range(min(10, len(ca_residues) - q_start)):
            pr = ca_residues[q_start + i][0]
            ur = uniprot_start_1based + i
            pa = pdb_seq[q_start + i] if (q_start + i) < len(pdb_seq) else "?"
            ua = uniprot_seq[t_start + i] if (t_start + i) < len(uniprot_seq) else "?"
            match = "OK" if pa == ua else "MISMATCH"
            print(f"  PDB {pr:4d} ({pa}) → UniProt {ur:4d} ({ua})  {match}")

        results[slug] = {
            "uniprot_accession": accession,
            "gene": gene,
            "offset": int(offset),
            "pdb_res_start": int(pdb_resSeq_start),
            "uniprot_res_start": int(uniprot_start_1based),
            "uniprot_res_end": int(uniprot_end_1based),
            "alignment_score": float(score),
            "coverage_pct": float(round(coverage, 1)),
        }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(str(OUT_JSON), "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nSaved: {OUT_JSON}")

    print("\n\nSUMMARY TABLE")
    print(f"{'Slug':<15} {'Gene':<8} {'Offset':>8} {'PDB_start':>10} "
          f"{'UP_start':>10} {'UP_end':>8} {'Cov%':>7} {'Score':>8}")
    print("-" * 80)
    for slug, r in results.items():
        print(
            f"{slug:<15} {r['gene']:<8} {r['offset']:>8} "
            f"{r['pdb_res_start']:>10} {r['uniprot_res_start']:>10} "
            f"{r['uniprot_res_end']:>8} {r['coverage_pct']:>6.1f}% "
            f"{r['alignment_score']:>8.0f}"
        )


if __name__ == "__main__":
    main()
