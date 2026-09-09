# SAMD9(L) Conservation and Mutational Landscape

An interactive figure comparing the amino-acid sequences of **SAMD9** and
**SAMD9L**, showing per-residue ConSurf conservation, domain boundaries, and
patient/population variants.

**Live site:** https://pmitev93.github.io/SAMD9_L_Map/

## Using the figure

- **Switch between Map and Table** with the buttons in the top-left. Map is the
  aligned-sequence figure; Table is every variant as a sortable list — click a
  column header to sort by it (click again to reverse), scroll and the header
  row stays frozen. Table also has its own **search box** (matches anything in
  the row), a **protein filter** (All / SAMD9 / SAMD9L), a **domain filter**,
  and a **conservation filter** (1-9, ConSurf score) — the last two are both
  multi-select checkbox panels (pick any combination of domains/linkers, or
  any combination of scores), independent of the category toggles below, all
  "everything" by default.
- **Export the table** with the CSV / Excel buttons above it — exports exactly
  the rows currently visible (search + protein/domain filter + category
  toggles + sort order all apply), not the whole dataset regardless of what
  you're looking at.
- **Toggle variant categories** with the box in the top-right (GoF / LoF / gnomAD
  / Somatic / NoF / Other / gnomAD missense). All gnomAD categories are off by
  default. These toggles filter **both** views — Table rows carry the same
  category as Map ticks, so hiding a category hides it everywhere.
- **Hover an amino-acid box** to see its exact position (e.g. `K133`).
- **Click a variant** to open a card with the source paper (linked), PMID,
  gnomAD status, phenotype, and method of functional assessment. gnomAD
  status is **live** — see below. The Table's gnomAD columns are also live,
  and include the population **homozygote count** (highlighted green when
  >0) alongside allele frequency.

## Files

| File | Purpose |
|------|---------|
| `index.html` | The figure itself (self-contained: CSS, fonts, legends and the domain scheme are all embedded). |
| `data_variants.js` | Every variant, one object each — base fields (protein, residue, label, effect) plus its popup details (paper, phenotype, method, gnomad) inline on the same object. |
| `data_papers.js` | Paper repository — each paper stored once, referenced by a key. |
| `data_overrides.js` | Manual label tweaks (nudge / line length). |
| `Conservation_Mutational_Landscape_Both.xlsx` | Master variant spreadsheet. |
| `tools/` | Build scripts + the renderer/CSS that get embedded into `index.html`. |

The three `data_*.js` files are **plain text you edit by hand**. After editing,
just **refresh the page** — no build step needed.

## How to…

### Add a variant
Add an entry to `data_variants.js` — one object holds everything for that
variant, both what renders it and what its popup shows:
```js
{
"protein": "SAMD9L",
"residue": 123,
"label": "K123R",
"effect": "GoF",
"origin": "smith_2020",
"paper": "smith_2020",
"phenotype": "GoF",
"method": "EdU assay"
},
```
`effect` sets the colour on the map: `GoF` (red), `LoF` (blue), `gnomAD`
(green), `Somatic` (black), `NoF` (light blue), anything else falls into
`Other` (gray) — see `deriveCategory()` in `tools/variant_renderer.js`.
`paper`/`phenotype`/`method` are optional — omit any of them and that popup
row just doesn't show, same as before. A variant with none of them shows
*"In progress"* in the popup — **except gnomAD status**, see below.
*Remember a comma after every `}` except the last one.*

### Add the paper info shown on click
Add the paper once in `data_papers.js`:
```js
"smith_2020": { "title": "Full paper title", "pmid": "12345678",
                "url": "https://pubmed.ncbi.nlm.nih.gov/12345678/" },
```
then reference it by key on the variant itself (`"paper": "smith_2020"`, as
in the example above) — no separate lookup file to keep in sync.

### gnomAD status — live, no manual entry needed
The page fetches gnomAD's variant data itself (2 requests on load, one per
gene, via gnomAD's public GraphQL API) and matches it to each variant by
protein position, so the popup's **gnomAD** row is always current — no field
to fill in. This covers missense (`R986C`), stop-gain (`W1507X`), frameshift
(`D1580VfsX2`), and single-residue deletion (`R1281del`) labels. Compound /
in-trans labels (e.g. `"R986C, T233N (in trans)"`) aren't a single
gnomAD-queryable variant, so they're skipped rather than guessed at.

You can still set a manual `"gnomad"` field on the variant itself in
`data_variants.js` — it's kept as a fallback shown only while the live fetch
is in flight, if it fails, or for labels live data can't speak to (compound
ones). No need to add it otherwise.

### gnomAD missense — auto-annotated, not just status-checked
Beyond checking gnomAD status on variants you've curated, the page also scans
gnomAD's live data for **missense** variants at residues with **no curated
annotation at all** (no GoF/LoF/Somatic/NoF/Other entry there) and adds them
as their own category — dark green, tick + label (e.g. `R986C`) just like
GoF/LoF — under the **gnomAD missense (unannotated)** toggle, off by default.
A residue that already has any curated variant is left as-is; gnomAD is never
used to second-guess an existing annotation. Truncating gnomAD variants
(stop-gain, frameshift, deletion) stay in the separate, tick-only **gnomAD**
category as before. This needs no manual upkeep — it's recomputed from the
same live fetch described above every time the page loads.

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
The Excel master only carries protein/residue/label/effect/origin — it knows
nothing about paper/phenotype/method. Before overwriting, this script reads
whatever's already in `data_variants.js` and re-attaches those detail fields
onto matching (protein, label) rows in the freshly-generated output, so a
rebuild doesn't wipe out hand-curated annotations. It prints a warning naming
any (protein, label) that had detail fields but no longer exists in the Excel
— those are genuinely lost (nothing to re-attach to) and need re-adding by
hand. Rows you've hand-added directly to `data_variants.js` that were never
in the Excel at all are dropped by a rebuild regardless (same as before) —
add those to the Excel master first, or re-add the whole row afterward.

### Deploy an update
Commit and push (e.g. with GitHub Desktop). The live site updates within ~1
minute. Hard-refresh (Cmd/Ctrl+Shift+R) if you don't see changes immediately.

### Change a domain boundary
Domain residue ranges live in ONE place, `DOMAIN_RANGES` near the top of
`tools/variant_renderer.js` — both the Map's colored outlines and the Table's
Domain column/filter/export read from it, so editing a boundary there updates
everywhere at once (then `python3 tools/embed_variants.py` to rebuild).

## Tuning knobs

Visual settings live in the `:root` "EASY TUNING KNOBS" block near the top of
`index.html` — e.g. `--vlabel-size`, `--vline-gap`, `--toggle-font`,
`--domain-height`, and the Walker-motif letter colours.
