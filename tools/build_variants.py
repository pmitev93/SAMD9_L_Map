#!/usr/bin/env python3
"""Convert the variant master spreadsheet into data_variants.js.

Reads Sheet1 of Conservation_Mutational_Landscape_Both.xlsx. Each data row is one
residue; the row's ID column is the residue number and the Mutation cell lists the
variant(s) at that residue. Distinct variants are separated by ';' and align with
';'-separated Effect/Origin; alternative notations within one variant use ','.

Output: data_variants.js — a flat list of
  {protein, residue, label, effect, origin, [group], [paper, gnomad, phenotype, method]}
(category is NOT stored — the renderer derives it from `effect` at load time,
see deriveCategory() in tools/variant_renderer.js) plus a small stats summary
printed to stdout (and any parse warnings).

paper/gnomad/phenotype/method are the per-variant popup details (this file
used to be split into data_variants.js + data_details.js; they were merged
into one file per variant). Since this script writes fresh rows straight from
the Excel master (which only carries protein/residue/label/effect/origin),
BEFORE overwriting it reads whatever data_variants.js already has on disk and
re-attaches any of those four detail fields onto matching (protein, label)
entries in the new output — so re-running this script does NOT wipe out
hand-curated paper/phenotype/method annotations the way it would if those
fields were only carried by the old Excel-derived row shape.

WARNING: this still OVERWRITES data_variants.js from the Excel master for the
base fields (protein/residue/label/effect/origin/group). If you have hand-
added a variant directly to data_variants.js that isn't in the Excel (e.g.
your own unpublished findings), running this WILL DISCARD that row entirely
(detail-preservation only helps for rows that still exist after the rebuild).
Add those to the Excel master first, or re-add the row afterward — the
detail fields you'd already set on it will NOT survive if the row itself is
dropped, only if its (protein, label) key still exists post-rebuild.

Run from the project folder:  python3 tools/build_variants.py
"""
import json, re, zipfile, pathlib, sys
from xml.etree import ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
XLSX = ROOT / "Conservation_Mutational_Landscape_Both.xlsx"
OUT  = ROOT / "data_variants.js"   # loaded via <script src> (works on double-click)
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# Column indices (1-based) in Sheet1, from the row-4 header.
COLS = {
    "SAMD9":  {"effect": 3, "mutation": 4, "id": 5, "origin": 2},
    "SAMD9L": {"id": 11, "mutation": 12, "effect": 13, "origin": 14},
}
FIRST_DATA_ROW = 5

# Mirrors deriveCategory() in tools/variant_renderer.js — kept here ONLY for the
# stdout stats summary below; it is not written into data_variants.js.
KNOWN_CATEGORIES = {"GoF", "LoF", "gnomAD", "Somatic", "NoF"}
def category(effect: str) -> str:
    first = re.split(r"[,/]", (effect or "").strip())[0].strip()
    return first if first in KNOWN_CATEGORIES else "Other"

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

DETAIL_FIELDS = ["paper", "gnomad", "phenotype", "method"]
BASE_FIELDS = ["protein", "residue", "label", "effect", "origin", "group"]

def load_existing_details():
    """Reads whatever data_variants.js already has on disk and returns
    {(protein, label): {detail-field: value}} for every entry that carries
    paper/gnomad/phenotype/method, so a rebuild can re-attach them to the
    freshly-generated rows below. "Last entry wins" per (protein, label),
    matching how the page itself resolves duplicates at runtime (DETAILS is
    built by iterating DATA in order, so a later duplicate's fields are what
    actually show in the popup)."""
    if not OUT.exists():
        return {}
    text = OUT.read_text(encoding="utf-8").strip()
    prefix = "window.VARIANTS ="
    if not text.startswith(prefix) or not text.endswith(";"):
        return {}
    try:
        existing = json.loads(text[len(prefix):-1].strip())
    except json.JSONDecodeError:
        return {}
    out = {}
    for v in existing:
        detail = {f_: v[f_] for f_ in DETAIL_FIELDS if v.get(f_) not in (None, "")}
        if detail:
            out[(v.get("protein"), v.get("label"))] = detail
    return out

def serialize(variants):
    """One field per line, no indentation — matches the hand-edited style in
    data_variants.js (base fields first, then any detail fields), so this
    script's output stays a small diff against manual edits."""
    lines = ["window.VARIANTS =", "["]
    for i, v in enumerate(variants):
        ordered = {f_: v[f_] for f_ in BASE_FIELDS if f_ in v}
        for f_ in DETAIL_FIELDS:
            if f_ in v:
                ordered[f_] = v[f_]
        lines.append("{")
        items = list(ordered.items())
        for j, (k, val) in enumerate(items):
            comma = "," if j < len(items) - 1 else ""
            lines.append(json.dumps(k) + ": " + json.dumps(val, ensure_ascii=False) + comma)
        lines.append("}" + ("," if i < len(variants) - 1 else ""))
    lines.append("];")
    return "\n".join(lines) + "\n"

def main():
    preserved = load_existing_details()
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
                                 "label": label, "effect": eff, "origin": ori})
    variants.sort(key=lambda v: (v["protein"], v["residue"], v["label"]))

    reattached = 0
    for v in variants:
        detail = preserved.get((v["protein"], v["label"]))
        if detail:
            v.update(detail)
            reattached += 1

    OUT.write_text(serialize(variants), encoding="utf-8")

    # ---- stats ----
    from collections import Counter
    by_prot = Counter(v["protein"] for v in variants)
    by_cat = Counter(category(v["effect"]) for v in variants)
    by_eff = Counter(v["effect"] for v in variants)
    print(f"Wrote {OUT.name}: {len(variants)} variants")
    print("  by protein:", dict(by_prot))
    print("  by category:", dict(by_cat))
    print("  by effect:", dict(by_eff))
    print(f"  re-attached paper/gnomad/phenotype/method to {reattached} rows "
          f"that already had them ({len(preserved)} distinct (protein,label) "
          f"keys carried detail fields before this rebuild)")
    dropped = set(preserved) - {(v["protein"], v["label"]) for v in variants}
    if dropped:
        print(f"\n  WARNING: {len(dropped)} (protein,label) keys had detail "
              f"fields but no longer exist in the Excel master — their "
              f"paper/gnomad/phenotype/method were DISCARDED:")
        for protein, label in sorted(dropped):
            print(f"   - {protein}:{label}")
    if warnings:
        print(f"\n  {len(warnings)} warnings (review):")
        for w in warnings[:25]:
            print("   -", w)
        if len(warnings) > 25:
            print(f"   ... +{len(warnings)-25} more")

if __name__ == "__main__":
    main()
