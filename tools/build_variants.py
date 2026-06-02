#!/usr/bin/env python3
"""Convert the variant master spreadsheet into variants.json.

Reads Sheet1 of Conservation_Mutational_Landscape_Both.xlsx. Each data row is one
residue; the row's ID column is the residue number and the Mutation cell lists the
variant(s) at that residue. Distinct variants are separated by ';' and align with
';'-separated Effect/Origin; alternative notations within one variant use ','.

Output: variants.json — a flat list of
  {protein, residue, label, effect, origin, category}
plus a small stats summary printed to stdout (and any parse warnings).

Run from the project folder:  python3 tools/build_variants.py
"""
import json, re, zipfile, pathlib, sys
from xml.etree import ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
XLSX = ROOT / "Conservation_Mutational_Landscape_Both.xlsx"
OUT  = ROOT / "variants.json"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# Column indices (1-based) in Sheet1, from the row-4 header.
COLS = {
    "SAMD9":  {"effect": 3, "mutation": 4, "id": 5, "origin": 2},
    "SAMD9L": {"id": 11, "mutation": 12, "effect": 13, "origin": 14},
}
FIRST_DATA_ROW = 5

def category(effect: str) -> str:
    # use the leading token (some cells join two effects with ',' or '/')
    first = re.split(r"[,/]", (effect or "").strip())[0].strip()
    if first == "GoF": return "GoF"
    if first == "LoF": return "LoF"
    if first == "gnomAD": return "gnomAD"
    return "Other"   # Somatic, NoF, blank, ...

def load_sheet1():
    z = ET.fromstring((zipfile.ZipFile(XLSX)).read("xl/sharedStrings.xml"))
    ss = ["".join(t.text or "" for t in si.iter(NS + "t")) for si in z]
    root = ET.fromstring(zipfile.ZipFile(XLSX).read("xl/worksheets/sheet1.xml"))
    rows = {}
    for c in root.iter(NS + "c"):
        m = re.match(r"([A-Z]+)(\d+)", c.get("r"))
        col = 0
        for ch in m.group(1):
            col = col * 26 + (ord(ch) - 64)
        row = int(m.group(2))
        v = c.find(NS + "v")
        if v is None:
            continue
        val = ss[int(v.text)] if c.get("t") == "s" else v.text
        rows.setdefault(row, {})[col] = val
    return rows

def split_aligned(mutation, effect, origin):
    """Split a cell into aligned (label, effect, origin) variant entries."""
    muts = [m.strip() for m in mutation.split(";")]
    effs = [e.strip() for e in (effect or "").split(";")]
    oris = [o.strip() for o in (origin or "").split(";")]
    out = []
    for i, label in enumerate(muts):
        if not label:
            continue
        eff = effs[i] if i < len(effs) else (effs[0] if effs else "")
        ori = oris[i] if i < len(oris) else (oris[0] if oris else "")
        out.append((label, eff, ori))
    return out, (len(muts) != len(effs))

def main():
    rows = load_sheet1()
    variants, warnings = [], []
    for r in sorted(rows):
        if r < FIRST_DATA_ROW:
            continue
        cells = rows[r]
        for protein, cmap in COLS.items():
            mut = cells.get(cmap["mutation"])
            if not mut or not str(mut).strip():
                continue
            rid = cells.get(cmap["id"])
            if rid is None:
                warnings.append(f"row {r} {protein}: mutation but no ID")
                continue
            residue = int(float(rid))
            entries, mism = split_aligned(str(mut), str(cells.get(cmap["effect"]) or ""),
                                          str(cells.get(cmap["origin"]) or ""))
            if mism:
                warnings.append(f"row {r} {protein} res {residue}: "
                                f"#mut!=#effect -> {mut!r} / {cells.get(cmap['effect'])!r}")
            for label, eff, ori in entries:
                variants.append({"protein": protein, "residue": residue,
                                 "label": label, "effect": eff, "origin": ori,
                                 "category": category(eff)})
    variants.sort(key=lambda v: (v["protein"], v["residue"], v["label"]))
    OUT.write_text(json.dumps(variants, ensure_ascii=False, indent=0), encoding="utf-8")

    # ---- stats ----
    from collections import Counter
    by_prot = Counter(v["protein"] for v in variants)
    by_cat = Counter(v["category"] for v in variants)
    by_eff = Counter(v["effect"] for v in variants)
    print(f"Wrote {OUT.name}: {len(variants)} variants")
    print("  by protein:", dict(by_prot))
    print("  by category:", dict(by_cat))
    print("  by effect:", dict(by_eff))
    if warnings:
        print(f"\n  {len(warnings)} warnings (review):")
        for w in warnings[:25]:
            print("   -", w)
        if len(warnings) > 25:
            print(f"   ... +{len(warnings)-25} more")

if __name__ == "__main__":
    main()
