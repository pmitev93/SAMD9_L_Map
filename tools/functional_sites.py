#!/usr/bin/env python3
"""For every curated variant, distance (in the AlphaFold model) to the
nearest known functional site — grounded in the REAL cryo-EM structures,
not guessed. Feeds the Table's "Nearest functional site" column.

Sites, both from Case Western's SAMD9 cryo-EM release (bioRxiv 2026.02.02.703423,
PDB 9ZJR/9ZJU — see project notes for how these were confirmed public):
  - nucleotide pocket:  SAMD9 residues (9ZJR, chain A) with any atom within
    SITE_RADIUS_A of the bound ATP/MG — this is the actual ligand density the
    paper reports ("ATP-bound... nucleotide-binding oligomerization domain"),
    not a guessed pocket.
  - dimer interface:    SAMD9 residues (9ZJU, chain A) with any atom within
    SITE_RADIUS_A of chain B — the real symmetric-dimer contact surface.
SAMD9L has no public cryo-EM structure yet (still embargoed — see project
notes), so its site residues are inferred by carrying the SAMD9 site
residues across data_residue_map.js's sequence correspondence. That's a
homology transfer, not a second independent structural observation — flagged
in the output and worth re-deriving once SAMD9L's own structure clears.

Distances themselves are measured on EACH protein's own AlphaFold model
(structures/SAMD9_AF.pdb / SAMD9L_AF_aligned.pdb) — the same coordinates the
3D panel already uses — not the cryo-EM structures directly, so a site
residue's identity comes from the real structure but its 3D position (and
every variant's) comes from the model already in use throughout the page.
For the ordered/well-resolved regions the two should track closely; this is
a documented approximation, not hidden.

Run from the project folder (after align_structures.py, which this depends
on for data_residue_map.js and SAMD9L_AF_aligned.pdb):
  python3 tools/functional_sites.py
"""
import json
import pathlib
import re
import numpy as np
from Bio.PDB import PDBParser

ROOT = pathlib.Path(__file__).resolve().parent.parent
CRYOEM_MONOMER = ROOT / "structures" / "cryoem_ref" / "9ZJR.pdb"   # SAMD9, ATP+Mg bound
CRYOEM_DIMER   = ROOT / "structures" / "cryoem_ref" / "9ZJU.pdb"   # SAMD9, symmetric dimer
AF_SAMD9       = ROOT / "structures" / "SAMD9_AF.pdb"
AF_SAMD9L      = ROOT / "structures" / "SAMD9L_AF_aligned.pdb"
RESIDUE_MAP_JS = ROOT / "data_residue_map.js"
VARIANTS_JS    = ROOT / "data_variants.js"
OUT_JS         = ROOT / "data_site_distances.js"

SITE_RADIUS_A = 5.0

SITE_LABELS = {"pocket": "nucleotide pocket", "interface": "dimer interface"}


def load_residue_ca(path, chain_id=None):
    """{resnum: np.array([x,y,z])} for a chain's CA atoms."""
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    chain = structure[0][chain_id] if chain_id else next(structure[0].get_chains())
    return {res.id[1]: res["CA"].coord for res in chain if "CA" in res}


def hetatm_coords(path, resnames, chain_id="A"):
    """All atom coordinates for the given HETATM residue name(s) on one chain."""
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    chain = structure[0][chain_id]
    coords = []
    for res in chain:
        if res.id[0].strip() and res.get_resname() in resnames:
            coords.extend(atom.coord for atom in res)
    return np.array(coords)


def residues_near(ca_map, target_coords, radius):
    """Residue numbers whose CA is within `radius` of ANY point in target_coords.
    CA-to-any-atom is a looser (larger) site than atom-to-atom would give —
    deliberately generous, since this defines "which residues count as the
    site" for a residue-level (not atom-level) distance metric downstream."""
    hits = set()
    for resnum, ca in ca_map.items():
        d = np.linalg.norm(target_coords - ca, axis=1)
        if d.min() <= radius:
            hits.add(resnum)
    return hits


def load_residue_map():
    text = RESIDUE_MAP_JS.read_text(encoding="utf-8")
    text = text[len("window.RESIDUE_MAP =\n"):].rstrip().rstrip(";")
    return json.loads(text)


def load_variants():
    text = VARIANTS_JS.read_text(encoding="utf-8").strip()
    prefix = "window.VARIANTS ="
    return json.loads(text[len(prefix):-1].strip())


def main():
    samd9_ca = load_residue_ca(AF_SAMD9)
    samd9l_ca = load_residue_ca(AF_SAMD9L)

    pocket_target = hetatm_coords(CRYOEM_MONOMER, {"ATP", "MG"}, "A")
    dimer_a_ca = load_residue_ca(CRYOEM_DIMER, "A")
    dimer_structure = PDBParser(QUIET=True).get_structure("d", str(CRYOEM_DIMER))
    interface_target = np.array([atom.coord for res in dimer_structure[0]["B"] for atom in res])

    samd9_pocket = residues_near(load_residue_ca(CRYOEM_MONOMER, "A"), pocket_target, SITE_RADIUS_A)
    samd9_interface = residues_near(dimer_a_ca, interface_target, SITE_RADIUS_A)
    print(f"SAMD9 (from real cryo-EM structures): {len(samd9_pocket)} nucleotide-pocket residues, "
          f"{len(samd9_interface)} dimer-interface residues")

    rmap = load_residue_map()
    fwd = rmap["SAMD9->SAMD9L"]
    samd9l_pocket = {fwd[str(r)] for r in samd9_pocket if str(r) in fwd}
    samd9l_interface = {fwd[str(r)] for r in samd9_interface if str(r) in fwd}
    print(f"SAMD9L (homology-transferred, no public structure yet): "
          f"{len(samd9l_pocket)} nucleotide-pocket, {len(samd9l_interface)} dimer-interface")

    sites = {
        "SAMD9":  {"pocket": (samd9_pocket, samd9_ca),  "interface": (samd9_interface, samd9_ca)},
        "SAMD9L": {"pocket": (samd9l_pocket, samd9l_ca), "interface": (samd9l_interface, samd9l_ca)},
    }

    variants = load_variants()
    seen = set()
    out = {}
    skipped = 0
    for v in variants:
        protein, residue = v["protein"], v["residue"]
        key = f"{protein}:{residue}"
        if key in seen:
            continue
        seen.add(key)
        ca_map = samd9_ca if protein == "SAMD9" else samd9l_ca
        if residue not in ca_map:
            skipped += 1
            continue
        my_ca = ca_map[residue]
        best_site, best_dist = None, None
        for site_key, (site_residues, site_ca_map) in sites[protein].items():
            for r in site_residues:
                if r == residue:
                    d = 0.0
                else:
                    d = float(np.linalg.norm(site_ca_map[r] - my_ca))
                if best_dist is None or d < best_dist:
                    best_dist, best_site = d, site_key
        if best_dist is not None:
            out[key] = {"distance": round(best_dist, 1), "site": SITE_LABELS[best_site]}

    js = ("window.SITE_DISTANCES =\n" + json.dumps(out, indent=0, sort_keys=True) + ";\n")
    OUT_JS.write_text(js, encoding="utf-8")
    print(f"Wrote {OUT_JS.relative_to(ROOT)}: {len(out)} residues "
          f"({len(seen)} distinct curated protein:residue keys, {skipped} skipped — "
          f"not found in that protein's AlphaFold model, e.g. a bad residue number)")


if __name__ == "__main__":
    main()
