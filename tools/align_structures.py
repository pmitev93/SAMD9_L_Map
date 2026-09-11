#!/usr/bin/env python3
"""Superimpose the SAMD9L AlphaFold model onto the SAMD9 one.

The 3D viewer (variant_renderer.js) shows one structure at a time and,
switching between them, now KEEPS the current camera (see renderResidue's
center()-based path) — so once the two structures share a coordinate frame,
switching proteins naturally lands on "the same view" of the corresponding
fold, without any extra viewer-side sync code. This script is what makes
them share that frame: it's a one-time (well, one-time-per-model-update)
precompute, not something the page redoes on every load.

Method: a global pairwise sequence alignment between the two canonical
sequences (BLOSUM62, affine gaps) gives candidate residue-residue
correspondences; their CA atoms seed a Kabsch/SVD superposition (Biopython's
SVDSuperimposer); one round of outlier pruning (drop pairs whose post-fit CA
distance sits far from the rest, then refit) keeps divergent loops from
skewing the core-fold fit — the two proteins are ~60% identical, not
identical, so a handful of pairs will always disagree, but the fold's rigid
core shouldn't. The resulting rotation+translation is applied to EVERY atom
of the SAMD9L model (not just CA) and written out as a new PDB; SAMD9's own
file is untouched (it's the reference frame).

Run from the project folder:  python3 tools/align_structures.py
"""
import json
import pathlib
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.Polypeptide import protein_letters_3to1_extended as AA3TO1
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF_PDB    = ROOT / "structures" / "SAMD9_AF.pdb"
MOVE_PDB   = ROOT / "structures" / "SAMD9L_AF.pdb"
OUT_PDB    = ROOT / "structures" / "SAMD9L_AF_aligned.pdb"
OUT_MAP_JS = ROOT / "data_residue_map.js"

OUTLIER_CUTOFF_A = 5.0   # post-fit CA distance beyond which a pair is dropped before refitting


def load_residues(path):
    """Returns (residue_numbers, one_letter_seq, ca_coords) in chain order."""
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    chain = next(structure[0].get_chains())
    resnums, seq, coords = [], [], []
    for res in chain:
        if "CA" not in res:
            continue
        resnums.append(res.id[1])
        seq.append(AA3TO1.get(res.get_resname(), "X"))
        coords.append(res["CA"].coord)
    return resnums, "".join(seq), np.array(coords, dtype=float)


def align_pairs(seq_ref, seq_move):
    """Global BLOSUM62 alignment -> list of (ref_index, move_index) for
    columns where both sides are a real (non-gap) residue."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(seq_ref, seq_move)[0]
    ref_idx, move_idx = aln.indices
    pairs = [(r, m) for r, m in zip(ref_idx, move_idx) if r != -1 and m != -1]
    return pairs


def fit(ref_coords, move_coords):
    sup = SVDSuperimposer()
    sup.set(ref_coords, move_coords)
    sup.run()
    rot, tran = sup.get_rotran()
    rms = sup.get_rms()
    return rot, tran, rms


def main():
    ref_resnums,  ref_seq,  ref_ca  = load_residues(REF_PDB)
    move_resnums, move_seq, move_ca = load_residues(MOVE_PDB)
    print(f"SAMD9  (reference): {len(ref_seq)} residues")
    print(f"SAMD9L (moving):    {len(move_seq)} residues")

    pairs = align_pairs(ref_seq, move_seq)
    print(f"Sequence alignment: {len(pairs)} aligned (non-gap) column pairs")

    ref_pts  = np.array([ref_ca[r]  for r, m in pairs])
    move_pts = np.array([move_ca[m] for r, m in pairs])

    rot, tran, rms = fit(ref_pts, move_pts)
    print(f"Initial fit: {len(pairs)} pairs, RMSD {rms:.2f} A")

    # One round of outlier rejection: refit without pairs whose CA lands far
    # from the reference after the initial (whole-alignment) superposition —
    # these are almost always divergent loops/termini, not the core fold.
    fitted = np.dot(move_pts, rot) + tran
    dist = np.linalg.norm(fitted - ref_pts, axis=1)
    keep = dist <= OUTLIER_CUTOFF_A
    dropped = int((~keep).sum())
    if dropped and keep.sum() >= 20:
        ref_pts2, move_pts2 = ref_pts[keep], move_pts[keep]
        rot, tran, rms = fit(ref_pts2, move_pts2)
        print(f"After dropping {dropped} outlier pairs (> {OUTLIER_CUTOFF_A} A): "
              f"{keep.sum()} pairs, RMSD {rms:.2f} A")
    else:
        print("No outlier-pruning refit needed/possible.")

    # Apply the final rotation+translation to EVERY atom (not just CA) of
    # the moving structure, and write it out.
    structure = PDBParser(QUIET=True).get_structure("SAMD9L", str(MOVE_PDB))
    for atom in structure.get_atoms():
        atom.transform(rot, tran)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(OUT_PDB))
    print(f"Wrote {OUT_PDB.relative_to(ROOT)}")

    # Residue-number correspondence for the "Compare to SAMD9(L)" feature
    # (variant_renderer.js) — the FULL sequence alignment (pairs, before the
    # outlier-pruning above), not just the pruned core used for the rigid
    # fit: pruning was about which pairs should influence the RIGID-BODY FIT,
    # not about which residues have a sequence correspondent at all. A loop
    # residue can be a perfectly good "analogous position" for comparison
    # even if its CA didn't superpose tightly.
    fwd, rev = {}, {}
    for r, m in pairs:
        fwd[ref_resnums[r]] = move_resnums[m]
        rev[move_resnums[m]] = ref_resnums[r]
    js = ("window.RESIDUE_MAP =\n" +
          json.dumps({"SAMD9->SAMD9L": fwd, "SAMD9L->SAMD9": rev}, indent=0) + ";\n")
    OUT_MAP_JS.write_text(js, encoding="utf-8")
    print(f"Wrote {OUT_MAP_JS.relative_to(ROOT)}: {len(fwd)} SAMD9->SAMD9L, "
          f"{len(rev)} SAMD9L->SAMD9 residue correspondences")


if __name__ == "__main__":
    main()
