# SAMD9(L) Conservation and Mutational Landscape

An interactive figure comparing the amino-acid sequences of **SAMD9** and
**SAMD9L**, showing per-residue ConSurf conservation, domain boundaries, and
patient/population variants.

**Live site:** https://pmitev93.github.io/SAMD9_L_Map/

## Using the figure

- **Toggle variant categories** with the box in the top-right (GoF / LoF / gnomAD
  / Other). gnomAD is off by default.
- **Hover an amino-acid box** to see its exact position (e.g. `K133`).
- **Click a variant** to open a card with the source paper (linked), PMID,
  gnomAD status, phenotype, and method of functional assessment. gnomAD
  status is **live** — see below.

## Files

| File | Purpose |
|------|---------|
| `index.html` | The figure itself (self-contained: CSS, fonts, legends and the domain scheme are all embedded). |
| `data_variants.js` | Every variant: protein, residue, label, effect, category. |
| `data_papers.js` | Paper repository — each paper stored once, referenced by a key. |
| `data_details.js` | Per-variant popup details, keyed `"PROTEIN:label"`. |
| `data_overrides.js` | Manual label tweaks (nudge / line length). |
| `Conservation_Mutational_Landscape_Both.xlsx` | Master variant spreadsheet. |
| `tools/` | Build scripts + the renderer/CSS that get embedded into `index.html`. |

The four `data_*.js` files are **plain text you edit by hand**. After editing,
just **refresh the page** — no build step needed.

## How to…

### Add a variant
Add an entry to `data_variants.js`:
```js
{"protein":"SAMD9L","residue":123,"label":"K123R","effect":"GoF","category":"GoF"},
```
`category` sets the colour: `GoF` (red), `LoF` (blue), `gnomAD` (green),
`Other` (black). *Remember a comma after every `}` except the last one.*

### Add the paper info shown on click
1. Add the paper once in `data_papers.js`:
```js
"smith_2020": { "title": "Full paper title", "pmid": "12345678",
                "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/" },
```
2. Link the variant to it in `data_details.js`:
```js
"SAMD9L:K123R": { "phenotype": "GoF", "method": "EdU assay", "paper": "smith_2020" },
```
Any field you omit is simply not shown. Variants without an entry (or without
a paper) show *"In progress"* — **except gnomAD status**, see below.

### gnomAD status — live, no manual entry needed
The page fetches gnomAD's variant data itself (2 requests on load, one per
gene, via gnomAD's public GraphQL API) and matches it to each variant by
protein position, so the popup's **gnomAD** row is always current — no field
to fill in. This covers missense (`R986C`), stop-gain (`W1507X`), frameshift
(`D1580VfsX2`), and single-residue deletion (`R1281del`) labels. Compound /
in-trans labels (e.g. `"R986C, T233N (in trans)"`) aren't a single
gnomAD-queryable variant, so they're skipped rather than guessed at.

You can still set a manual `"gnomad"` field in `data_details.js` — it's kept
as a fallback shown only while the live fetch is in flight, if it fails, or
for labels live data can't speak to (compound ones). No need to add it
otherwise.

If gnomAD ever ships a new default release, `tools/variant_renderer.js` has
two constants (`GNOMAD_DATASET`, `GNOMAD_VERSION_LABEL`) that need a manual
bump — check the version shown in the page title at
[gnomad.broadinstitute.org](https://gnomad.broadinstitute.org). The allele
frequencies themselves update automatically regardless; only that label
needs occasional attention.

### Nudge a label (overlap / crowding)
In `data_overrides.js`, keyed by label:
```js
"V1276I": { "dx": -8, "dlen": 10 }
```
- `dx` — move the label left(−)/right(+)
- `dlen` — make the line longer(+)/shorter(−); the label follows
- `lane` — snap to a stacking level (`0` = shortest)

### Re-build variants from the spreadsheet
Only needed if you edit the Excel master rather than `data_variants.js`
directly (requires Python 3):
```
python3 tools/build_variants.py
```

### Deploy an update
Commit and push (e.g. with GitHub Desktop). The live site updates within ~1
minute. Hard-refresh (Cmd/Ctrl+Shift+R) if you don't see changes immediately.

## Tuning knobs

Visual settings live in the `:root` "EASY TUNING KNOBS" block near the top of
`index.html` — e.g. `--vlabel-size`, `--vline-gap`, `--toggle-font`,
`--domain-height`, and the Walker-motif letter colours.
