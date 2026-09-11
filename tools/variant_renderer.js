/* Data-driven variant renderer for the SAMD9/SAMD9L map.
   Reads embedded variant-data (+ optional variant-overrides), removes the old
   hand-placed markers, and draws colored ticks + lane-stacked labels above
   SAMD9L / below SAMD9. Also builds the show/hide category toggles. */
(function () {
  "use strict";

  // ---- Toggle/legend text: edit the "legend" strings below to rename a category.
  var CFG = {
    GoF:       { color: "#FF0000", label: true,  on: true,  legend: "Gain-of-function" },
    LoF:       { color: "#1f56bc", label: true,  on: true,  legend: "Loss-of-function" },
    gnomAD:    { color: "#73d73c", label: false, on: false, legend: "gnomAD (truncating)" },
    Somatic:   { color: "#000000", label: true,  on: false, legend: "Somatic" },
    NoF:       { color: "#8ECAE6", label: true,  on: false, legend: "NoF (no gain-of-function effect)" },
    Other:     { color: "#888888", label: true,  on: false, legend: "Other" },
    // Missense variants seen in gnomAD's live data at a residue we haven't
    // otherwise annotated (see injectGnomadMissenseVariants below) — labeled
    // like GoF/LoF (unlike the tick-only truncating gnomAD category above),
    // in a darker green so the two stay visually distinct.
    gnomADmis: { color: "#1B6B32", label: true,  on: false, legend: "gnomAD missense (unannotated)" }
  };
  // category is DERIVED from `effect` (not stored) — one source of truth. The
  // leading token before a "," or "/" decides the bucket; anything unrecognized
  // (e.g. "SIRT2") falls into the Other catch-all so it still renders.
  function deriveCategory(effect) {
    var first = (effect || "").split(/[,/]/)[0].trim();
    if (CFG[first]) return first;
    return "Other";
  }

  // ---- Live gnomAD data — fetched from gnomAD's own public API on every page
  // load, so no one has to manually look up and re-enter allele frequencies.
  // Queries the WHOLE gene once per protein (not one request per variant),
  // filtered to the canonical/MANE-Select transcript, and matches each
  // variant's protein-level HGVS (e.g. "p.Arg986Cys") back to our own label
  // format ("R986C") so it can be looked up directly by data-mutation.
  //
  // GNOMAD_DATASET / GNOMAD_VERSION_LABEL: gnomAD's API has no queryable
  // "current version" field, so the display label is a plain string here —
  // if gnomAD ships a new default release, bump GNOMAD_DATASET (e.g. to
  // "gnomad_r5") and GNOMAD_VERSION_LABEL to match (check the version shown
  // in the page title at gnomad.broadinstitute.org). The allele-frequency
  // numbers themselves stay live/automatic regardless; only this label needs
  // an occasional manual nudge.
  var GNOMAD_DATASET = "gnomad_r4";
  var GNOMAD_VERSION_LABEL = "gnomAD v4.1.1";
  var GNOMAD_GENES = [
    { protein: "SAMD9",  symbol: "SAMD9",  transcript: "ENST00000379958" },  // MANE Select
    { protein: "SAMD9L", symbol: "SAMD9L", transcript: "ENST00000318238" }   // MANE Select
  ];
  var GNOMAD_LIVE = { ready: false, failed: false, byLabel: { SAMD9: {}, SAMD9L: {} },
                       missense: { SAMD9: [], SAMD9L: [] } };
  var AA3TO1 = {
    Ala: "A", Arg: "R", Asn: "N", Asp: "D", Cys: "C", Gln: "Q", Glu: "E", Gly: "G",
    His: "H", Ile: "I", Leu: "L", Lys: "K", Met: "M", Phe: "F", Pro: "P", Ser: "S",
    Thr: "T", Trp: "W", Tyr: "Y", Val: "V", Ter: "X"
  };
  // Converts gnomAD's protein HGVS into our own label convention, covering the
  // shapes that actually occur in this dataset: missense, stop-gain,
  // frameshift, single-residue deletion. Anything else (splice, extension,
  // multi-residue indels...) returns null and is simply not matched.
  function hgvspToLabel(hgvsp) {
    if (!hgvsp) return null;
    var m;
    if ((m = hgvsp.match(/^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})fsTer(\d+)$/))) {
      var ref1 = AA3TO1[m[1]], alt1 = AA3TO1[m[3]];
      return (ref1 && alt1) ? ref1 + m[2] + alt1 + "fsX" + m[4] : null;
    }
    if ((m = hgvsp.match(/^p\.([A-Za-z]{3})(\d+)del$/))) {
      var ref2 = AA3TO1[m[1]];
      return ref2 ? ref2 + m[2] + "del" : null;
    }
    if ((m = hgvsp.match(/^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$/))) {
      var ref3 = AA3TO1[m[1]], alt3 = AA3TO1[m[3]];
      return (ref3 && alt3) ? ref3 + m[2] + alt3 : null;
    }
    return null;
  }
  // A variant label counts as "askable" of gnomAD only if it's a single,
  // simple protein change in one of the shapes above — compound/in-trans
  // labels (e.g. "R986C, T233N (in trans)") aren't a single gnomAD variant,
  // so we never claim to know whether gnomAD "has" them.
  function isSimpleLabel(label) {
    return /^[A-Z]\d+([A-Z]|del|fsX\d+|[A-Z]fsX\d+)$/.test(label);
  }
  function fetchGnomadGene(symbol) {
    var query = "query($sym:String!){ gene(gene_symbol:$sym, reference_genome: GRCh38) { " +
                "variants(dataset: " + GNOMAD_DATASET + ") { hgvsp transcript_id " +
                "exome{ac an af homozygote_count} genome{ac an af homozygote_count} } } }";
    return fetch("https://gnomad.broadinstitute.org/api", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ query: query, variables: { sym: symbol } })
    }).then(function (r) { return r.json(); });
  }
  // A label counts as strict missense (not stop-gain, frameshift, or
  // deletion, OR synonymous) when it's Ref+Pos+Alt, Alt is a real amino acid
  // (not the "X" this project uses for a stop — hgvspToLabel funnels
  // stop-gain through this same Ref+Pos+Alt shape, Ter -> "X"), and Alt
  // actually differs from Ref (Ref===Alt is a synonymous/silent change —
  // same codon position mutated, same amino acid — not a missense variant).
  function isMissenseLabel(label) {
    var m = label.match(/^([A-Z])(\d+)([A-Z])$/);
    return !!m && m[3] !== "X" && m[1] !== m[3];
  }
  function loadGnomadLive() {
    Promise.all(GNOMAD_GENES.map(function (g) {
      return fetchGnomadGene(g.symbol).then(function (json) {
        var variants = (json.data && json.data.gene && json.data.gene.variants) || [];
        variants.forEach(function (v) {
          if (v.transcript_id !== g.transcript) return;
          var label = hgvspToLabel(v.hgvsp);
          if (!label) return;
          var freq = v.exome || v.genome;
          if (!freq) return;
          GNOMAD_LIVE.byLabel[g.protein][label] = { ac: freq.ac, an: freq.an, af: freq.af, hom: freq.homozygote_count };
          if (freq.ac > 0 && isMissenseLabel(label)) {
            var residue = +label.match(/^[A-Z](\d+)/)[1];
            GNOMAD_LIVE.missense[g.protein].push({ label: label, residue: residue });
          }
        });
      });
    })).then(function () {
      GNOMAD_LIVE.ready = true;
      injectGnomadMissenseVariants();
    }).catch(function () {
      GNOMAD_LIVE.failed = true;
    });
  }
  // Adds a tick+label (like GoF/LoF) for any gnomAD missense variant at a
  // residue we have NOT already annotated with something else (any existing
  // curated variant at that residue — any category — wins; we never add a
  // duplicate or second opinion at an already-annotated position). Re-runs
  // the renderer once so the new entries actually appear, and adds the
  // category's checkbox to the already-built toggle box.
  var gnomadMissenseInjected = false;
  function injectGnomadMissenseVariants() {
    if (gnomadMissenseInjected) return;   // guard: loadGnomadLive only ever resolves once
    gnomadMissenseInjected = true;
    ["SAMD9", "SAMD9L"].forEach(function (protein) {
      var annotatedResidues = {};
      DATA.forEach(function (v) { if (v.protein === protein) annotatedResidues[v.residue] = true; });
      GNOMAD_LIVE.missense[protein].forEach(function (m) {
        if (annotatedResidues[m.residue]) return;   // already annotated -> leave as-is
        DATA.push({
          protein: protein, residue: m.residue, label: m.label,
          effect: "gnomADmis", category: "gnomADmis", origin: "gnomAD (live)"
        });
        annotatedResidues[m.residue] = true;   // gnomAD can list >1 missense per residue; keep only the first
      });
    });
    addToggleFor("gnomADmis");
    run();
    buildTableView();
  }
  // Formats the gnomAD popup row. Live data wins whenever we have it (a real
  // match, or a confirmed absence from the full canonical-transcript variant
  // list); falls back to the manually-curated "gnomad" field (on the variant
  // itself, in data_variants.js) only when live data can't speak to this
  // specific label (still loading, fetch failed, or the label is a compound
  // variant gnomAD can't be asked about).
  function gnomadRowFor(protein, label, staticGnomad) {
    if (GNOMAD_LIVE.ready && isSimpleLabel(label)) {
      var hit = GNOMAD_LIVE.byLabel[protein] && GNOMAD_LIVE.byLabel[protein][label];
      if (hit && hit.ac > 0) {
        return "Present — allele frequency " + hit.af.toExponential(2) +
               " (" + hit.ac + "/" + hit.an + ") — " + GNOMAD_VERSION_LABEL;
      }
      return "Not present — " + GNOMAD_VERSION_LABEL;
    }
    if (staticGnomad != null) {
      var g = (typeof staticGnomad === "object")
        ? (staticGnomad.present ? "Yes" + (staticGnomad.maf ? " (MAF " + staticGnomad.maf + ")" : "") : "No")
        : staticGnomad;
      return g;
    }
    if (GNOMAD_LIVE.failed) return null;
    return isSimpleLabel(label) ? "Loading live gnomAD data…" : null;
  }
  // Geometry. LANE_H/LABEL_H track the label size so bigger labels still fit.
  // TICK = the (fixed) line length for a lane-0 variant; +LANE_H per stack level.
  var TICK = 26, LANE_H = 15, LABEL_H = 15, PAD = 5;
  // white space between a line's bottom and the AA box (the --vline-gap knob)
  function whiteGap() {
    var v = parseFloat(getComputedStyle(document.documentElement)
              .getPropertyValue("--vline-gap"));
    return isNaN(v) ? 5 : v;
  }

  function estW(s) { return s.length * 7.0 + 8; }   // ~ label pixel width

  function hasDigit(el) { return el && /\d/.test(el.textContent || ""); }

  // Find the number-ruler row next to a sequence row, skipping blank spacer rows
  // (some blocks have an empty <tr> between the sequence and its ruler). dir = +1
  // looks below (SAMD9), -1 looks above (SAMD9L). Returns the ruler row or null.
  function findRuler(row, dir) {
    var n = dir > 0 ? row.nextElementSibling : row.previousElementSibling;
    for (var s = 0; n && s < 4; s++) {
      var t = (n.textContent || "").trim();
      if (/\d/.test(t)) return n;            // the ruler
      if (t !== "" && t !== ".") return null; // hit real content -> no ruler here
      n = dir > 0 ? n.nextElementSibling : n.previousElementSibling;
    }
    return null;
  }

  // ---- Domain boundaries + color utilities — single source of truth shared
  // by the Map's gradient outlines (run(), below) AND the Table's Domain
  // column/filter/export. N-to-C order matters: it's how linker/terminus
  // names between consecutive domains get generated. Boundaries were
  // extracted directly from the live, previously-correct page (not
  // re-derived/guessed), residue by residue, so they match the original
  // domains exactly.
  var DOMAIN_RANGES = [
    { name: "SAM",           key: "sam",     SAMD9L: [10, 86],     SAMD9: [10, 86] },
    { name: "AlbA",          key: "alba",    SAMD9L: [157, 389],   SAMD9: [157, 384] },
    { name: "SIR2",          key: "sir2",    SAMD9L: [394, 625],   SAMD9: [390, 622] },
    { name: "P-loop NTPase", key: "ploop",   SAMD9L: [712, 1017],  SAMD9: [708, 1013] },
    { name: "TPR",           key: "tpr",     SAMD9L: [1027, 1187], SAMD9: [1023, 1187] },
    { name: "Helical",       key: "helical", SAMD9L: [1193, 1497], SAMD9: [1193, 1502] },
    { name: "OB-fold",       key: "obfold",  SAMD9L: [1511, 1579], SAMD9: [1516, 1584] }
  ];
  function hexToRgb(hex) {
    hex = hex.replace("#", "");
    if (hex.length === 3) hex = hex.split("").map(function (c) { return c + c; }).join("");
    return [parseInt(hex.slice(0, 2), 16), parseInt(hex.slice(2, 4), 16), parseInt(hex.slice(4, 6), 16)];
  }
  // Named CSS colors (yellow/blue/orange/red/purple) don't parse via hexToRgb
  // directly — resolve anything through the browser first.
  function toHex(cssColor) {
    var probe = document.createElement("div");
    probe.style.color = cssColor;
    document.body.appendChild(probe);
    var rgb = getComputedStyle(probe).color.match(/\d+/g).map(Number);
    document.body.removeChild(probe);
    return "#" + rgb.map(function (v) { return v.toString(16).padStart(2, "0"); }).join("");
  }
  function lerpColor(hexA, hexB, t) {
    var a = hexToRgb(hexA), b = hexToRgb(hexB);
    var rgb = a.map(function (v, i) { return Math.round(v + (b[i] - v) * t); });
    return "rgb(" + rgb.join(",") + ")";
  }
  // Blend a color toward white by `amount` (0-1) to get its "light" partner.
  function lighten(hex, amount) {
    var rgb = hexToRgb(hex).map(function (v) { return Math.round(v + (255 - v) * amount); });
    return "#" + rgb.map(function (v) { return v.toString(16).padStart(2, "0"); }).join("");
  }
  function knob(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
  }
  // Oscillating wave, not a straight sweep: color alternates base -> light ->
  // base -> light... every WAVELENGTH residues, all the way across the domain
  // (continues smoothly across row-wraps since position is domain-relative).
  var WAVELENGTH = 14; // residues per full dark->light->dark cycle — tune here
  function waveColor(baseHex, lightHex, i) {
    var phase = (i % WAVELENGTH) / WAVELENGTH;
    var osc = (1 - Math.cos(2 * Math.PI * phase)) / 2;   // 0..1..0 smoothly
    return lerpColor(baseHex, lightHex, osc);
  }
  // Resolved once (cheap, but no reason to redo it per row): a flat
  // representative color per domain — for the Table's badges and exports,
  // NOT the Map's animated dark/light wave. Reads the same --domain-*-dark
  // knobs run() uses, so retuning those stays in sync everywhere. Linkers and
  // the N-/C-termini get one shared neutral gray (they're not a domain).
  var DOMAIN_COLOR_CACHE = null;
  function domainColors() {
    if (DOMAIN_COLOR_CACHE) return DOMAIN_COLOR_CACHE;
    DOMAIN_COLOR_CACHE = {
      sam:     toHex(knob("--domain-sam-dark", "yellow")),
      alba:    toHex(knob("--domain-alba-dark", "blue")),
      sir2:    toHex(knob("--domain-sir2-dark", "orange")),
      ploop:   toHex(knob("--domain-ploop-dark", "red")),
      tpr:     toHex(knob("--domain-tpr-dark", "#16D715")),
      obfold:  toHex(knob("--domain-obfold-dark", "purple")),
      helical: toHex(knob("--domain-helical-dark", "cyan")),
      linker:  "#9AA0A6", nterm: "#9AA0A6", cterm: "#9AA0A6", none: "#9AA0A6"
    };
    return DOMAIN_COLOR_CACHE;
  }
  // Which domain a residue falls in — or, for a residue between two domains
  // (or before the first / after the last), a synthesized "X–Y linker" /
  // "N-terminus" / "C-terminus" name so every residue gets SOME label.
  function domainInfoFor(protein, residue) {
    for (var i = 0; i < DOMAIN_RANGES.length; i++) {
      var range = DOMAIN_RANGES[i][protein];
      if (range && residue >= range[0] && residue <= range[1]) {
        return { name: DOMAIN_RANGES[i].name, key: DOMAIN_RANGES[i].key };
      }
    }
    var prev = null, next = null;
    for (var j = 0; j < DOMAIN_RANGES.length; j++) {
      var r = DOMAIN_RANGES[j][protein];
      if (!r) continue;
      if (residue < r[0] && !next) next = DOMAIN_RANGES[j];
      if (residue > r[1]) prev = DOMAIN_RANGES[j];
    }
    if (prev && next) return { name: prev.name + "–" + next.name + " linker", key: "linker" };
    if (!prev && next) return { name: "N-terminus", key: "nterm" };
    if (prev && !next) return { name: "C-terminus", key: "cterm" };
    return { name: "—", key: "none" };
  }
  // Canonical N-to-C ordering for the Table's domain filter dropdown —
  // domains interleaved with the linker names domainInfoFor() actually
  // produces, so every option in the list can actually match a row.
  var DOMAIN_FILTER_ORDER = (function () {
    var order = ["N-terminus"];
    DOMAIN_RANGES.forEach(function (d, i) {
      order.push(d.name);
      if (i < DOMAIN_RANGES.length - 1) order.push(d.name + "–" + DOMAIN_RANGES[i + 1].name + " linker");
    });
    order.push("C-terminus");
    return order;
  })();

  // ConSurf conservation score (1-9) per residue — read straight off each
  // sequence cell's own "ScoreN" CSS class (run() sets RESIDUE_CELLS once,
  // below; the residue->cell mapping is static so one capture covers every
  // later Table build). Alignment-gap dashes carry "Score0", which is not a
  // real conservation score — no residue label ever lands on a gap cell, so
  // that value is never actually returned here.
  var RESIDUE_CELLS = null;
  function conservationFor(protein, residue) {
    var cell = RESIDUE_CELLS && RESIDUE_CELLS[protein] && RESIDUE_CELLS[protein][residue];
    if (!cell) return null;
    var m = cell.className.match(/Score(\d+)/);
    return m ? +m[1] : null;
  }

  // Distance (Å, in that protein's AlphaFold model) from a curated variant's
  // residue to the nearest known functional site — precomputed offline from
  // the REAL SAMD9 cryo-EM structures (tools/functional_sites.py), not
  // guessed; see that script's own header for exactly how, and for the
  // SAMD9L caveat (homology-transferred, no public SAMD9L structure yet).
  // { distance, site } or null (e.g. an alignment-gap residue with no
  // AlphaFold coordinate at all).
  function siteFor(protein, residue) {
    return (window.SITE_DISTANCES && window.SITE_DISTANCES[protein + ":" + residue]) || null;
  }

  // ---- data: read from the data_*.js globals (loaded via <script src>, which
  //      works even on a double-clicked file). Edit a data_*.js + refresh = live.
  var DATA = [], OV = {}, PAPERS = {}, DETAILS = {};
  function ovFor(v) { return OV[v.label] || OV[v.protein + ":" + v.label] || null; }

  function init() {
    DATA    = window.VARIANTS || [];
    DATA.forEach(function (v) { v.category = deriveCategory(v.effect); });
    OV      = window.VARIANT_OVERRIDES || {};
    PAPERS  = window.PAPERS || {};
    // DETAILS used to be a separate data_details.js file keyed "PROTEIN:label".
    // paper/gnomad/phenotype/method now live directly on each DATA entry (one
    // file, one object per variant) — DETAILS is just that same lookup,
    // rebuilt from DATA itself, so openPopup()/buildTableView() (which read
    // d.paper / d.phenotype / d.method / d.gnomad off whatever this returns)
    // don't need to change at all.
    DETAILS = {};
    DATA.forEach(function (v) { DETAILS[v.protein + ":" + v.label] = v; });
    run();
    buildViewSwitcher();
    buildTableView();
    setupPopup();
  }

  function run() {
    var table = document.getElementById("samd9-table");
    if (!table) return;
    // Skip entirely while the Map is hidden (Table view active): every
    // measurement below is a getBoundingClientRect, which returns all-zero
    // rects for a display:none subtree, and the debounced window "resize"
    // listener (schedule() at the bottom of this file) can absolutely fire
    // while hidden — the Browser pane detaching/reattaching alone does it.
    // Without this guard that corrupts every .vlayer's computed height into
    // garbage (confirmed: 450-560px instead of the normal ~30-60px) that
    // then persists after switching back to Map, since nothing else
    // recomputes it. setView() forces one fresh run() on Map re-activation
    // to guarantee correct geometry even if the window was resized while
    // Table was showing.
    if (table.offsetParent === null) return;

    // 1. drop legacy markers + any previous render (idempotent on resize)
    table.querySelectorAll(".vertical-line").forEach(function (n) { n.remove(); });
    table.querySelectorAll("tr.vmark").forEach(function (n) { n.remove(); });
    // 2. drop the now-empty old marker rows (top/bottom-align rows with no text).
    //    Number rulers keep their digits and survive; sequence rows aren't
    //    top/bottom-align, so they're untouched.
    table.querySelectorAll("tr.top-align, tr.bottom-align").forEach(function (r) {
      if (!(r.textContent || "").trim()) r.remove();
    });

    // 3. insert a reserved-height container OUTSIDE the number ruler, so the
    //    ruler stays glued to the sequence (above SAMD9L / below SAMD9).
    var seqRows = [];
    Array.prototype.forEach.call(table.querySelectorAll("tr"), function (row) {
      var first = row.querySelector("td");
      if (!first) return;
      var t = (first.textContent || "").trim();
      if (t !== "SAMD9" && t !== "SAMD9L") return;
      var side = (t === "SAMD9L") ? "top" : "bottom";
      var tr = document.createElement("tr");
      tr.className = "vmark vmark-" + side;
      var td = document.createElement("td");
      td.colSpan = 60;
      var div = document.createElement("div");
      div.className = "vlayer";
      td.appendChild(div); tr.appendChild(td);
      if (side === "top") {
        var anchor = findRuler(row, -1) || row;          // ruler above SAMD9L
        anchor.parentNode.insertBefore(tr, anchor);
      } else {
        var a2 = findRuler(row, 1) || row;               // ruler below SAMD9
        a2.parentNode.insertBefore(tr, a2.nextSibling);
      }
      seqRows.push({ protein: t, row: row, side: side, container: div });
    });

    // 4. residue -> cell map (cumulative per protein, document order). Also
    // track EVERY sequence cell (letters AND alignment-gap "-" dashes) per
    // protein, in document order, with an index — lets domain outlines below
    // span a residue range without leaving a gap-dash cell unstyled.
    var maps = { SAMD9: {}, SAMD9L: {} };
    var allCells = { SAMD9: [], SAMD9L: [] };
    var cellIndex = new Map();
    var cellRow = new Map();
    var counters = { SAMD9: 0, SAMD9L: 0 };
    seqRows.forEach(function (sr) {
      var cells = sr.row.children;
      for (var i = 1; i < cells.length; i++) {
        var txt = (cells[i].textContent || "").trim();
        var isLetter = /^[A-Z]$/.test(txt), isGap = txt === "-";
        cellRow.set(cells[i], sr);
        // allCells only holds genuine sequence positions (a real residue
        // letter, or an alignment-gap dash). Some rows have a truly blank
        // leading spacer cell (empty text, class="Score11") between the
        // protein-name label and the first residue — including that here
        // would let a domain outline spill onto it and color a cell that
        // isn't part of the sequence at all.
        if (isLetter || isGap) {
          cellIndex.set(cells[i], allCells[sr.protein].length);
          allCells[sr.protein].push(cells[i]);
        }
        if (isLetter) {
          counters[sr.protein]++;
          maps[sr.protein][counters[sr.protein]] = cells[i];
          cells[i].dataset.pos = txt + counters[sr.protein];   // e.g. "K133"
          cells[i].dataset.side = sr.side;                     // top / bottom
          cells[i].dataset.protein = sr.protein;                // SAMD9 / SAMD9L
          cells[i].dataset.resnum = counters[sr.protein];       // numeric residue #, for the 3D-structure click (setupResidueClick, below)
        }
      }
    });
    // Residue -> cell never changes between runs (the underlying sequence
    // markup is static; only the variant overlay changes) — cache once for
    // the Table's Conservation column/filter/export (conservationFor(),
    // below) so it doesn't need its own copy of this scan.
    RESIDUE_CELLS = maps;

    // 4b. domain-boundary outlines (data-driven). All 7 domains render through
    // this one path — the older SAM/AlbA/SIR2/P-loop NTPase/TPR/OB-fold used to
    // be hand-coded per-row .outlined-row-all-* CSS classes (still present in
    // the stylesheet, now inert/superseded — inline styles here always win).
    // Migrating them here is what makes a TRUE horizontal gradient possible:
    // each cell's color is computed from its own fractional position across
    // the whole domain and applied identically to top AND bottom, so the two
    // lines can never mismatch. Boundaries below were extracted directly from
    // the live, previously-correct page (not re-derived/guessed), residue by
    // residue, so they match the original domains exactly. (The color
    // utilities and the boundaries themselves now live at module scope —
    // shared with the Table's Domain column/filter/export — see
    // DOMAIN_RANGES / domainColors() above.)
    var PAL = domainColors();
    var DOMAIN_OUTLINES = [];
    DOMAIN_RANGES.forEach(function (d) {
      ["SAMD9L", "SAMD9"].forEach(function (protein) {
        var range = d[protein];
        if (range) DOMAIN_OUTLINES.push({ protein: protein, start: range[0], end: range[1], pal: PAL[d.key] });
      });
    });
    var DOMAIN_BORDER_W = 2; // px — thinner than the original 3px
    DOMAIN_OUTLINES.forEach(function (d) {
      var startCell = maps[d.protein] && maps[d.protein][d.start];
      var endCell   = maps[d.protein] && maps[d.protein][d.end];
      if (!startCell || !endCell) return;
      var idxStart = cellIndex.get(startCell), idxEnd = cellIndex.get(endCell);
      var cells = allCells[d.protein].slice(idxStart, idxEnd + 1);
      var base = d.pal, light = lighten(d.pal, 0.55);
      cells.forEach(function (cell, i) {
        var color = waveColor(base, light, i);   // same color for top+bottom -> lines always match
        cell.style.borderTop = DOMAIN_BORDER_W + "px solid " + color;
        cell.style.borderBottom = DOMAIN_BORDER_W + "px solid " + color;
        cell.style.borderLeft = (i === 0 ? DOMAIN_BORDER_W + "px" : "0px") + " solid " + color;
        cell.style.borderRight = (i === cells.length - 1 ? DOMAIN_BORDER_W + "px" : "0px") + " solid " + color;
        // (tried a black keyline outside the colored border via box-shadow —
        // dropped it: table row paint order means the row above/below often
        // paints over the shadow's peeking sliver, so it was invisible for
        // most of the domain and looked broken at the caps. Not worth the
        // fragility; removed.)
        // (also tried rounded "pill" caps at the true start/end — dropped:
        // read as trimmed/oval rather than a clean square box. Plain square
        // corners, matching the original domains.)
      });
    });

    // 5. group variants by block
    var groups = new Map(), missing = 0;
    DATA.forEach(function (v) {
      var cell = maps[v.protein] && maps[v.protein][v.residue];
      if (!cell) { missing++; return; }
      v._cell = cell;
      var sr = cellRow.get(cell);
      if (!groups.has(sr)) groups.set(sr, []);
      groups.get(sr).push(v);
    });

    // 6. lay out + draw
    groups.forEach(function (list, sr) {
      var base = sr.container.getBoundingClientRect().left;
      list.forEach(function (v) {
        var r = v._cell.getBoundingClientRect();
        v._x = r.left - base + r.width / 2;
      });
      list.sort(function (a, b) { return a._x - b._x || a.residue - b.residue; });

      list.forEach(function (v) {
        v._color = (CFG[v.category] || CFG.Other).color;
        var o = ovFor(v);
        if (o && o.color) v._color = o.color;   // per-variant color override
      });

      // Cluster same-residue variants sharing an explicit `group` id into one
      // visual unit: one tick, labels placed side-by-side on the same line
      // (comma-separated), each still independently clickable/colored. A
      // variant without `group` is its own 1-item cluster (unchanged behavior).
      var clusters = [], seenGroups = {};
      list.forEach(function (v) {
        if (v.group) {
          if (seenGroups[v.group]) return;
          seenGroups[v.group] = true;
          clusters.push(list.filter(function (o) { return o.group === v.group; }));
        } else {
          clusters.push([v]);
        }
      });

      var lanes = [], maxLane = -1;
      clusters.forEach(function (c) {
        var rep = c[0];
        // Only members whose OWN category shows a label (e.g. not gnomAD) take
        // part in the visible line's text/width. An unlabeled member (e.g. a
        // gnomAD co-mutation compounded with a labeled one) still gets its own
        // tick drawn — see drawCluster — it just contributes no text/width here,
        // so it can never push a sibling's label off-center or hide it via the
        // category toggle (each tick is tagged with its OWN member's category).
        var labeled = c.filter(function (m) { return (CFG[m.category] || CFG.Other).label; });
        labeled.forEach(function (m, i) {
          m._text = m.label + (i < labeled.length - 1 ? ", " : "");
          m._w = estW(m._text);
        });
        var totalW = labeled.reduce(function (s, m) { return s + m._w; }, 0);
        var l = -1;
        if (labeled.length) {
          var leftEdge = rep._x - totalW / 2; l = 0;
          while (l < lanes.length && lanes[l] > leftEdge - 3) l++;
          lanes[l] = rep._x + totalW / 2;
        }
        c.forEach(function (m) { m._lane = (labeled.indexOf(m) >= 0) ? l : -1; });
      });
      // apply manual overrides: dx = horizontal label nudge; lane = stack level;
      // dlen = extra line length in px (label moves out with the line).
      list.forEach(function (v) {
        var o = ovFor(v);
        v._dx = (o && o.dx) || 0;
        v._dlen = (o && o.dlen) || 0;
        if (o && o.lane != null && v._lane >= 0) v._lane = o.lane;
        if (v._lane > maxLane) maxLane = v._lane;
      });

      // distance from the container edge to the sequence letters. A number ruler
      // (when present) sits in this gap; it varies block to block, so we use it
      // only to anchor the whitespace — NOT the line length (kept fixed below).
      var cr = sr.container.getBoundingClientRect(), rr = sr.row.getBoundingClientRect();
      var gap = Math.max(0, Math.round(sr.side === "top" ? rr.top - cr.bottom
                                                         : cr.top - rr.bottom));
      sr._gap = gap;
      // reserve only the part of the tallest stack that rises above the ruler band
      var W = whiteGap();
      var maxLen = TICK + (maxLane >= 0 ? maxLane * LANE_H : 0);
      list.forEach(function (v) {                       // account for any dlen extensions
        if (v._lane >= 0) maxLen = Math.max(maxLen, TICK + v._lane * LANE_H + v._dlen);
      });
      var topAboveBox = W + maxLen + (maxLane >= 0 ? LABEL_H : 0) + PAD;
      sr.container.style.height = Math.max(0, topAboveBox - gap) + "px";
      clusters.forEach(function (c) { drawCluster(sr, c); });
    });

    buildToggles();
    setupPosTip(table);
    setupResidueClick(table);
    window.__variantInfo = { total: DATA.length, missing: missing, counts: counters };
    if (missing) console.warn("variant renderer: " + missing + " unmapped variants");
  }

  // Draws one cluster: one tick PER MEMBER (grouped members share a residue,
  // so same-category ticks simply coincide — looks like one line — while a
  // mixed-category group, e.g. a gnomAD co-mutation alongside a labeled one,
  // shows each color independently and each tick is only ever hidden by ITS
  // OWN category's toggle, never a sibling's). Plus one `.vlabel` per LABELED
  // member, laid out side-by-side on the same line so grouped variants (e.g.
  // "R1281K, R1281S, R1281del") read as one line while each label stays its
  // own independently clickable/colored element. A plain (ungrouped) variant
  // is just a 1-member cluster — identical to the old draw().
  function drawCluster(sr, c) {
    var top = (sr.side === "top");
    var gap = sr._gap || 0;
    var W = whiteGap();              // white space before the AA box
    var base = gap - W;              // line bottom sits W above the actual letter box

    c.forEach(function (m) {
      var len = TICK + (m._lane >= 0 ? m._lane * LANE_H : 0) + (m._dlen || 0);
      var tick = document.createElement("div");
      tick.className = "vtick";
      tick.style.left = (m._x - 1) + "px";
      tick.style.height = len + "px";
      tick.style.background = m._color;
      tick.style.color = m._color;   // for the hover glow (currentColor)
      tick.style[top ? "bottom" : "top"] = (-base) + "px";
      setData(tick, m);
      sr.container.appendChild(tick);
    });

    var labeled = c.filter(function (m) { return m._lane >= 0; });
    if (!labeled.length) return;

    var rep = labeled[0];
    var len = TICK + rep._lane * LANE_H + (rep._dlen || 0);
    var totalW = labeled.reduce(function (s, m) { return s + m._w; }, 0);
    var leftEdge = rep._x - totalW / 2 + (rep._dx || 0);
    var running = 0;
    labeled.forEach(function (m) {
      var cx = leftEdge + running + m._w / 2 + (m === rep ? 0 : (m._dx || 0));
      running += m._w;
      var lab = document.createElement("div");
      lab.className = "vlabel";
      lab.textContent = m._text;
      lab.style.left = cx + "px";
      lab.style.color = m._color;
      lab.style[top ? "bottom" : "top"] = (-base + len + 1) + "px";
      setData(lab, m);
      sr.container.appendChild(lab);
    });
  }

  function setData(el, v) {
    el.setAttribute("data-category", v.category);
    el.setAttribute("data-protein", v.protein);
    el.setAttribute("data-residue", v.residue);
    el.setAttribute("data-mutation", v.label);
    el.setAttribute("data-origin", v.origin || "");
    el.title = v.label + " — click for details";
  }

  // ---- click-a-variant popup (glassy card with the paper / gnomAD info) ----
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }
  function setupPopup() {
    if (window.__vpop) return;
    var pop = document.createElement("div");
    pop.id = "variant-popup"; pop.style.display = "none";
    document.body.appendChild(pop);
    window.__vpop = pop;
    document.addEventListener("click", function (e) {
      var m = e.target.closest && e.target.closest(".vtick,.vlabel");
      if (m) { e.stopPropagation(); openPopup(m, pop); }
      // A variant sphere in the 3D viewer (renderAllVariants, below) opens
      // its OWN popup via a 3Dmol shape callback that fires on the
      // underlying 'mouseup', one tick before the browser's own synthetic
      // 'click' on the same interaction reaches here — without this
      // exclusion that later click would immediately re-close what the
      // callback just opened. Any other click inside either viewer (empty
      // space, a plain atom, a drag) just leaves the popup alone rather
      // than closing it.
      else if (!(e.target.closest && (e.target.closest("#variant-popup") ||
                 e.target.closest("#sp-viewer, #sp-viewer-2")))) hide(pop);
    });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") hide(pop); });
    window.addEventListener("resize", function () { hide(pop); });
  }
  // keep the clicked label's popup-owner in front of same-line neighbors even
  // after the mouse leaves it (e.g. moves onto the popup itself)
  var activeMarker = null;
  function hide(pop) {
    pop.style.display = "none";
    if (activeMarker) { activeMarker.classList.remove("vlabel-active"); activeMarker = null; }
  }

  function openPopup(marker, pop) {
    if (activeMarker) activeMarker.classList.remove("vlabel-active");
    activeMarker = marker.classList.contains("vlabel") ? marker : null;
    if (activeMarker) activeMarker.classList.add("vlabel-active");
    var protein = marker.getAttribute("data-protein");
    var mut = marker.getAttribute("data-mutation");
    var cat = marker.getAttribute("data-category");
    var origin = marker.getAttribute("data-origin");
    var color = (CFG[cat] || CFG.Other).color;
    var d = DETAILS[protein + ":" + mut] || DETAILS[mut];

    var h = '<button class="vpop-close" aria-label="Close">&times;</button>';
    h += '<div class="vpop-head"><span class="vpop-chip" style="background:' + color +
         '"></span>' + esc(mut) + '<span class="vpop-prot">' + esc(protein) + '</span></div>';

    var anyRow = false;
    if (d) {
      var paper = (d.paper && PAPERS[d.paper]) || {};
      var url = paper.url || (paper.pmid ? "https://pubmed.ncbi.nlm.nih.gov/" + paper.pmid + "/" : null);
      if (paper.title) {
        h += '<div class="vpop-row"><b>Source:</b> ' +
             '<a class="vpop-title" ' + (url ? 'href="' + esc(url) + '" target="_blank" rel="noopener"' : "") +
             ' title="' + esc(paper.title) + '">' + esc(paper.title) + "</a></div>";
        anyRow = true;
      }
      if (paper.pmid) {
        h += '<div class="vpop-row vpop-dim"><b>PMID:</b> ' +
             '<a href="https://pubmed.ncbi.nlm.nih.gov/' + esc(paper.pmid) + '/" target="_blank" rel="noopener">' +
             esc(paper.pmid) + "</a></div>";
        anyRow = true;
      }
    }
    // gnomAD: independent of `d` — live data (or its "loading"/absence state)
    // shows for ANY variant, not just ones with a manually-curated entry.
    var gnomadText = gnomadRowFor(protein, mut, d && d.gnomad);
    if (gnomadText != null) {
      h += '<div class="vpop-row"><b>gnomAD:</b> ' + esc(gnomadText) + "</div>";
      anyRow = true;
    }
    if (d) {
      var pheno = d.phenotype || d.effect;
      if (pheno) { h += '<div class="vpop-row"><b>Phenotype:</b> ' + esc(pheno) + "</div>"; anyRow = true; }
      if (d.method) {
        h += '<div class="vpop-row"><b>Method of functional assessment:</b> ' + esc(d.method) + "</div>";
        anyRow = true;
      }
    }
    if (!anyRow) h += '<div class="vpop-row vpop-dim"><i>In progress</i></div>';
    pop.innerHTML = h;
    pop.style.display = "block";
    pop.querySelector(".vpop-close").onclick = function () { hide(pop); };
    positionPopup(pop, marker);
  }

  function positionPopup(pop, marker) {
    var r = marker.getBoundingClientRect();
    var pw = pop.offsetWidth, ph = pop.offsetHeight, M = 8;
    var left = window.scrollX + r.left + r.width / 2 - pw / 2;
    left = Math.max(window.scrollX + M, Math.min(left, window.scrollX + document.documentElement.clientWidth - pw - M));
    var below = r.bottom + ph + M < window.innerHeight;
    var top = window.scrollY + (below ? r.bottom + M : r.top - ph - M);
    pop.style.left = left + "px";
    pop.style.top = Math.max(window.scrollY + M, top) + "px";
  }

  // Hover an amino-acid box -> show its exact position (e.g. "K133"), above the
  // box for SAMD9L, below it for SAMD9. Plain black text, no background.
  function setupPosTip(table) {
    if (table.__posTip) return;
    table.__posTip = true;
    var tip = document.createElement("div");
    tip.id = "aa-postip";
    document.body.appendChild(tip);
    table.addEventListener("mouseover", function (e) {
      var cell = e.target.closest ? e.target.closest("td[data-pos]") : null;
      if (!cell) return;
      var r = cell.getBoundingClientRect();
      tip.textContent = cell.dataset.pos;
      tip.style.display = "block";
      tip.style.left = (window.scrollX + r.left + r.width / 2) + "px";
      var above = cell.dataset.side === "top";
      tip.style.top = (window.scrollY + (above ? r.top - 16 : r.bottom + 4)) + "px";
    });
    table.addEventListener("mouseout", function (e) {
      if (e.target.closest && e.target.closest("td[data-pos]")) tip.style.display = "none";
    });
  }

  // ---- 3D structure viewer (AlphaFold model, opens on any residue click) ----
  // A right-side panel that shows the clicked residue's position in that
  // protein's AlphaFold model (3Dmol.js). Deliberately lazy in both pieces:
  // 3Dmol.js itself (~600KB) and each protein's structure file (~1MB PDB,
  // local under structures/ — not fetched from AlphaFold DB at runtime, so
  // this doesn't depend on their CORS policy or uptime) only ever load the
  // first time they're actually needed, so a visitor who never clicks a
  // residue pays nothing for this feature. Step 1 of a planned series (see
  // project notes) — AlphaFold only for now, one structure at a time, no
  // superposition/neighbor-highlighting/export yet.
  // SAMD9L's file is a SUPERPOSED copy (tools/align_structures.py), not the
  // raw AlphaFold download — coordinates pre-rotated/translated into
  // SAMD9's own frame so the two share one coordinate system. SAMD9 is the
  // fixed reference and needs no transform. That's what lets the camera
  // (renderResidue's center()-based path) carry over an orientation across
  // a protein switch and actually land on the corresponding view, not just
  // the same zoom level. Re-run that script if either AlphaFold model updates.
  var STRUCTURE_SOURCES = {
    SAMD9:  { file: "structures/SAMD9_AF.pdb" },
    SAMD9L: { file: "structures/SAMD9L_AF_aligned.pdb" }
  };
  // Same 1-9 ConSurf palette as the 2D map's td.ScoreN cells (index.html's
  // inline <style>) — kept as a literal copy, not read off the DOM, since
  // there's no single shared cartoon-vs-td color source to read from.
  var CONSERVATION_COLORS = {
    1: "#107f84", 2: "#44afbf", 3: "#a6dde7", 4: "#d7eef2", 5: "#FFFFFF",
    6: "#fbecf4", 7: "#f9c9dd", 8: "#f07dab", 9: "#a12561"
  };
  // Deliberately NOT the Table's domainColors().linker/nterm/cterm/none grey
  // (#9AA0A6) — that one's tuned to read on the Table's white background;
  // this one only ever sits in the 3D viewer (black by default), so it gets
  // its own, lighter value instead of the two use-cases fighting over one.
  var STRUCTURE_NO_DOMAIN_GREY = "#D8DBDF";
  // 3-letter PDB residue name -> 1-letter code, for labeling whatever the
  // user clicks directly ON the structure (residueClickLabel, below) in the
  // same "K620" style the rest of the page already uses.
  var AA_3TO1 = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
    GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
    PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V"
  };
  var structureText = {};   // protein -> already-fetched PDB text (session cache)
  var mol3dReady = null;    // becomes a resolved Promise once 3Dmol.js has loaded
  var structureViewer = null;   // one persistent $3Dmol.GLViewer, reused across clicks
  var loadedProtein = null;     // which protein's model is currently addModel()'d into it
  var currentResnum = null;     // last-highlighted residue, for re-applying style on a mode/bg toggle
  var structureColorMode = "conservation";   // "conservation" | "domain"
  var structureBg = "white";                 // "black" | "white"
  var hasFramedView = false;    // true after the first-ever render — gates the one-time zoomTo
                                 // vs. every render after that (incl. switching protein) using
                                 // center(), so a zoom level the user set up survives a protein switch
  var isFullscreen = false;
  var panelWidthPct = 50;       // restored when exiting fullscreen

  // ---- Compare: second, synced-camera viewer showing the OTHER protein at
  // the analogous residue (data_residue_map.js, from tools/align_structures.py) ----
  var compareMode = false;
  var structureViewer2 = null;
  var loadedProtein2 = null;
  var currentResnum2 = null;

  function otherProtein(protein) { return protein === "SAMD9" ? "SAMD9L" : "SAMD9"; }
  // Residue-number correspondence is a byproduct of the SAME sequence
  // alignment align_structures.py runs to superpose the two structures —
  // exported alongside it as window.RESIDUE_MAP so this lookup is just a
  // dictionary read, no alignment logic duplicated in the browser. Returns
  // null for a position with no 1:1 correspondent (an indel column).
  function analogousResidue(protein, resnum) {
    var map = window.RESIDUE_MAP && window.RESIDUE_MAP[protein + "->" + otherProtein(protein)];
    var mapped = map && map[resnum];
    return mapped != null ? mapped : null;
  }

  function load3Dmol() {
    if (window.$3Dmol) return Promise.resolve();
    if (mol3dReady) return mol3dReady;
    mol3dReady = new Promise(function (resolve, reject) {
      var s = document.createElement("script");
      s.src = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js";
      s.onload = function () { resolve(); };
      s.onerror = function () { mol3dReady = null; reject(new Error("couldn't load 3Dmol.js")); };
      document.head.appendChild(s);
    });
    return mol3dReady;
  }

  function fetchStructure(protein) {
    if (structureText[protein]) return Promise.resolve(structureText[protein]);
    var src = STRUCTURE_SOURCES[protein];
    if (!src) return Promise.reject(new Error("no structure file for " + protein));
    return fetch(src.file).then(function (r) {
      if (!r.ok) throw new Error("HTTP " + r.status);
      return r.text();
    }).then(function (text) {
      structureText[protein] = text;
      return text;
    });
  }

  // One [data-value] segmented control per group ("color": conservation/
  // domain, "bg": black/white) — shared markup+listener, dispatched by
  // data-group below. Doesn't touch the camera; only setColorMode/
  // setStructureBackground do, and neither of those moves it either.
  function segHtml(group, options, active) {
    return '<div class="sp-seg" data-group="' + group + '">' +
      options.map(function (o) {
        return '<button type="button" class="sp-seg-btn' + (o.value === active ? " sp-seg-active" : "") +
               '" data-value="' + o.value + '">' + o.label + "</button>";
      }).join("") +
      "</div>";
  }

  function buildStructurePanel() {
    if (document.getElementById("structure-panel")) return;
    document.documentElement.style.setProperty("--sp-width", panelWidthPct + "%");
    var panel = document.createElement("div");
    panel.id = "structure-panel";
    panel.innerHTML =
      '<div class="sp-drag" id="sp-drag" title="Drag to resize"></div>' +
      '<div class="sp-head">' +
        '<span class="sp-title" id="sp-title">Structure</span>' +
        '<button type="button" class="sp-close" id="sp-close" aria-label="Close">&times;</button>' +
      "</div>" +
      '<div class="sp-controls">' +
        segHtml("color", [{ value: "conservation", label: "Conservation" }, { value: "domain", label: "Domain" }], structureColorMode) +
        segHtml("bg", [{ value: "black", label: "Black" }, { value: "white", label: "White" }], structureBg) +
        '<button type="button" class="sp-fs-btn" id="sp-neighbors-btn">Nearby Residues</button>' +
        '<button type="button" class="sp-fs-btn" id="sp-allvariants-btn">All Variants</button>' +
        '<button type="button" class="sp-fs-btn" id="sp-compare-btn">Compare</button>' +
        '<button type="button" class="sp-fs-btn" id="sp-fs-btn">Fullscreen</button>' +
      "</div>" +
      '<div class="sp-body">' +
        '<div class="sp-status" id="sp-status"></div>' +
        '<div class="sp-viewer-row" id="sp-viewer-row">' +
          '<div class="sp-viewer-pane" id="sp-viewer-pane">' +
            '<div id="sp-viewer"></div>' +
            '<div class="sp-pane-cap" id="sp-pane-cap-1"></div>' +
          "</div>" +
          '<div class="sp-viewer-pane" id="sp-viewer-pane-2" style="display:none;">' +
            '<div id="sp-viewer-2"></div>' +
            '<div class="sp-pane-cap" id="sp-pane-cap-2"></div>' +
          "</div>" +
        "</div>" +
        '<div class="sp-clickhint" id="sp-clickhint">Click an atom on the structure to identify it</div>' +
      "</div>" +
      '<div class="sp-neighbors" id="sp-neighbors" style="display:none;"></div>' +
      '<div class="sp-foot">AlphaFold model — predicted structure, not experimental.</div>';
    document.body.appendChild(panel);
    panel.querySelector("#sp-close").addEventListener("click", closeStructurePanel);
    panel.querySelector("#sp-fs-btn").addEventListener("click", toggleFullscreen);
    panel.querySelector("#sp-neighbors-btn").addEventListener("click", toggleNeighbors);
    panel.querySelector("#sp-allvariants-btn").addEventListener("click", toggleAllVariants);
    panel.querySelector("#sp-compare-btn").addEventListener("click", toggleCompare);
    panel.querySelector(".sp-controls").addEventListener("click", function (e) {
      var btn = e.target.closest && e.target.closest(".sp-seg-btn");
      if (!btn) return;
      var seg = btn.closest(".sp-seg");
      seg.querySelectorAll(".sp-seg-btn").forEach(function (b) { b.classList.toggle("sp-seg-active", b === btn); });
      var group = seg.getAttribute("data-group"), value = btn.getAttribute("data-value");
      if (group === "color") setColorMode(value);
      else if (group === "bg") setStructureBackground(value);
    });
    setupDrag(panel.querySelector("#sp-drag"));
    wireCompareSync();
  }

  function toggleFullscreen() {
    isFullscreen = !isFullscreen;
    var panel = document.getElementById("structure-panel");
    panel.classList.toggle("sp-fullscreen", isFullscreen);
    document.getElementById("sp-fs-btn").textContent = isFullscreen ? "Exit Fullscreen" : "Fullscreen";
    resizeStructureViewer();
  }

  // Throttled to one resize+render per animation frame — a drag fires many
  // mousemove events, and a WebGL canvas resize isn't free enough to do on
  // every single one of them.
  var resizePending = false;
  function resizeStructureViewer() {
    if (!structureViewer || resizePending) return;
    resizePending = true;
    requestAnimationFrame(function () {
      resizePending = false;
      structureViewer.resize();
      structureViewer.render();
      if (structureViewer2) { structureViewer2.resize(); structureViewer2.render(); }
    });
  }

  // Drag the panel's left edge to resize it against the map — updates the
  // shared --sp-width custom property that both the panel and the map's own
  // shrunk-width rules (variant_styles.css) read, so they stay in sync with
  // one write instead of coordinating two separate elements.
  function setupDrag(handle) {
    var dragging = false;
    handle.addEventListener("mousedown", function (e) {
      if (isFullscreen) return;
      dragging = true;
      document.body.style.cursor = "col-resize";
      document.body.style.userSelect = "none";
      e.preventDefault();
    });
    document.addEventListener("mousemove", function (e) {
      if (!dragging) return;
      var pct = 100 - (e.clientX / window.innerWidth * 100);
      pct = Math.max(20, Math.min(85, pct));
      panelWidthPct = pct;
      document.documentElement.style.setProperty("--sp-width", pct + "%");
      resizeStructureViewer();
    });
    document.addEventListener("mouseup", function () {
      if (!dragging) return;
      dragging = false;
      document.body.style.cursor = "";
      document.body.style.userSelect = "";
    });
  }

  // Re-style only — never touches the camera, so toggling color/background
  // mid-inspection can't undo a zoom/pan the user set up (see renderResidue).
  function setColorMode(mode) {
    structureColorMode = mode;
    if (!structureViewer || !loadedProtein) return;
    applyStructureStyle(loadedProtein);
    if (currentResnum != null) highlightResidue(currentResnum);
    if (showingNeighbors) renderNeighbors(); else structureViewer.render();
    if (compareMode && structureViewer2 && loadedProtein2) {
      applyStructureStyle2(loadedProtein2);
      if (currentResnum2 != null) highlightResidue2(currentResnum2);
      structureViewer2.render();
    }
  }
  function setStructureBackground(color) {
    structureBg = color;
    if (structureViewer) { structureViewer.setBackgroundColor(color); structureViewer.render(); }
    if (structureViewer2) { structureViewer2.setBackgroundColor(color); structureViewer2.render(); }
  }

  // ---- Nearby residues (5A) + polar contacts ----
  var NEIGHBOR_RADIUS = 5;      // Angstrom — "nearby" residues
  var POLAR_DISTANCE  = 3.5;    // Angstrom — heteroatom pair counted as a polar contact
  var POLAR_ELEMS     = { N: 1, O: 1 };
  var BACKBONE_ATOMS  = { N: 1, CA: 1, C: 1, O: 1 };
  var showingNeighbors = false;

  function dist3(a, b) {
    var dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  function toggleNeighbors() {
    showingNeighbors = !showingNeighbors;
    document.getElementById("sp-neighbors-btn").classList.toggle("sp-fs-btn-active", showingNeighbors);
    // Mutually exclusive with "All Variants" — both draw into the same
    // addShape/addLabel buckets on the SAME viewer, so showing both at once
    // would just be visual noise on top of each other, not a combined view.
    if (showingNeighbors && showingAllVariants) toggleAllVariants();
    if (showingNeighbors) renderNeighbors(); else clearNeighbors();
  }

  function clearNeighbors() {
    if (structureViewer && loadedProtein) {
      structureViewer.removeAllShapes();
      structureViewer.removeAllLabels();
      applyStructureStyle(loadedProtein);
      if (currentResnum != null) highlightResidue(currentResnum);
      structureViewer.render();
    }
    var text = document.getElementById("sp-neighbors");
    if (text) { text.innerHTML = ""; text.style.display = "none"; }
  }

  // ---- All curated variants on the structure at once, colored like the Map ----
  var showingAllVariants = false;

  function toggleAllVariants() {
    showingAllVariants = !showingAllVariants;
    document.getElementById("sp-allvariants-btn").classList.toggle("sp-fs-btn-active", showingAllVariants);
    if (showingAllVariants && showingNeighbors) toggleNeighbors();
    if (showingAllVariants) renderAllVariants(); else clearAllVariants();
  }

  function clearAllVariants() {
    if (structureViewer && loadedProtein) {
      structureViewer.removeAllShapes();
      applyStructureStyle(loadedProtein);
      if (currentResnum != null) highlightResidue(currentResnum);
      structureViewer.render();
    }
  }

  // A sphere per curated variant on the CURRENT protein (its own category
  // color — same CFG the Map/legend use), skipping whatever categories are
  // currently hidden via the category toggles so this stays consistent with
  // the Map/Table rather than a separate, disconnected filter state. Each
  // sphere is independently clickable — opens the SAME popup a Map tick
  // would (openVariantPopupAt, below), not a simplified stand-in.
  function renderAllVariants() {
    if (!structureViewer || !loadedProtein) return;
    var protein = loadedProtein;
    structureViewer.removeAllShapes();
    DATA.filter(function (v) {
      return v.protein === protein && !document.body.classList.contains("hide-cat-" + v.category);
    }).forEach(function (v) {
      var atom = structureViewer.selectedAtoms({ resi: v.residue, atom: "CA" })[0];
      if (!atom) return;
      var color = (CFG[v.category] || CFG.Other).color;
      structureViewer.addSphere({
        center: { x: atom.x, y: atom.y, z: atom.z }, radius: .9, color: color,
        clickable: true,
        callback: function (shape, viewer, event) {
          if (event && event.stopPropagation) event.stopPropagation();
          openVariantPopupAt(v, event ? event.clientX : window.innerWidth / 2, event ? event.clientY : window.innerHeight / 2);
        }
      });
    });
    structureViewer.render();
  }

  // Opens the exact same popup a Map tick/label does (openPopup, setupPopup
  // above), for a click that didn't originate on a real DOM marker — a
  // variant sphere in the 3D view. openPopup only needs an element with the
  // right data-* attributes and a getBoundingClientRect(); a zero-size,
  // invisible div positioned at the click's screen coordinates satisfies
  // that with no changes to openPopup itself, then is discarded — openPopup
  // reads it synchronously (position, attributes) and never needs it again.
  function openVariantPopupAt(variant, clientX, clientY) {
    var pop = window.__vpop;
    if (!pop) return;
    var anchor = document.createElement("div");
    anchor.style.cssText = "position:fixed; left:" + clientX + "px; top:" + clientY + "px; width:0; height:0;";
    anchor.setAttribute("data-protein", variant.protein);
    anchor.setAttribute("data-mutation", variant.label);
    anchor.setAttribute("data-category", variant.category);
    anchor.setAttribute("data-origin", variant.origin || "");
    document.body.appendChild(anchor);
    openPopup(anchor, pop);
    anchor.remove();
  }

  // Residues with any atom within NEIGHBOR_RADIUS of the selected residue,
  // plus polar contacts: a SIDE-CHAIN N/O atom on the selected residue
  // within POLAR_DISTANCE of an N/O atom on a DIFFERENT residue (backbone
  // atoms allowed on that other side — the "side chain" restriction in the
  // ask is about OUR residue's own side chain, not the partner's).
  // Finds every curated variant at a residue (there can be more than one —
  // e.g. R620Q and R620W both at 620) — used to decide whether a "nearby
  // residue" pill (renderNeighbors, below) is clickable.
  function variantsAt(protein, resnum) {
    return DATA.filter(function (v) { return v.protein === protein && v.residue === resnum; });
  }
  // One pill: plain if this residue has no curated variant, or — reusing
  // the SAME class + data-* attributes setData() puts on a Map tick/label —
  // a real ".vlabel" if it does, so the page's EXISTING document-level click
  // listener (setupPopup()) opens the normal variant popup for it with no
  // extra wiring here. (Multiple variants at one residue: only the first is
  // linked — a known simplification, rare in this dataset.)
  function neighborPillHtml(protein, resnum, resn, isSelf) {
    var score = conservationFor(protein, resnum);
    var color = (score != null && CONSERVATION_COLORS[score]) || "#bbb";
    var text = (AA_3TO1[resn] || resn) + resnum;
    var variant = variantsAt(protein, resnum)[0];
    var style = "--pill-c:" + color;
    if (isSelf) return '<span class="sp-pill sp-pill-self" style="' + style + '">' + esc(text) + "</span>";
    if (!variant) return '<span class="sp-pill" style="' + style + '">' + esc(text) + "</span>";
    return '<span class="sp-pill vlabel" style="' + style + '" data-protein="' + esc(variant.protein) +
      '" data-mutation="' + esc(variant.label) + '" data-category="' + esc(variant.category) +
      '" data-origin="' + esc(variant.origin || "") + '" title="Click to view ' + esc(variant.label) + '">' +
      esc(text) + "</span>";
  }

  function renderNeighbors() {
    if (!structureViewer || !loadedProtein || currentResnum == null) return;
    var protein = loadedProtein, resnum = currentResnum;
    var targetAtoms = structureViewer.selectedAtoms({ resi: resnum });
    if (!targetAtoms.length) return;
    var nearAtoms = structureViewer.selectedAtoms({ within: { distance: NEIGHBOR_RADIUS, sel: { resi: resnum } } })
      .filter(function (a) { return a.resi !== resnum; });

    // One representative atom per nearby residue (prefer CA — a sensible,
    // stable anchor point for that residue's in-scene label below).
    var nearbyReps = {};
    nearAtoms.forEach(function (a) {
      if (!nearbyReps[a.resi] || a.atom === "CA") nearbyReps[a.resi] = a;
    });
    var nearbyResnums = Object.keys(nearbyReps).map(Number).sort(function (a, b) { return a - b; });

    var targetPolar = targetAtoms.filter(function (a) { return !BACKBONE_ATOMS[a.atom] && POLAR_ELEMS[a.elem]; });
    var nearPolar = nearAtoms.filter(function (a) { return POLAR_ELEMS[a.elem]; });
    var polarResidues = {}, lines = [];
    targetPolar.forEach(function (ta) {
      nearPolar.forEach(function (na) {
        if (dist3(ta, na) <= POLAR_DISTANCE) {
          lines.push([ta, na]);
          polarResidues[na.resi] = na.resn;
        }
      });
    });
    var polarResnums = Object.keys(polarResidues).map(Number).sort(function (a, b) { return a - b; });

    applyStructureStyle(protein);
    highlightResidue(resnum);
    structureViewer.removeAllShapes();
    structureViewer.removeAllLabels();
    if (nearbyResnums.length) {
      structureViewer.addStyle({ resi: nearbyResnums }, { stick: { colorscheme: "Jmol", radius: .13 } });
    }
    // addLine's dashed mode is a real WebGL line — most browsers/GPUs clamp
    // gl.lineWidth to 1px regardless of the requested width, so it barely
    // shows. addCylinder({dashed:true}) draws actual thin dashed CYLINDERS
    // instead, genuinely thicker via `radius`, not subject to that clamp.
    lines.forEach(function (pair) {
      structureViewer.addCylinder({
        start: { x: pair[0].x, y: pair[0].y, z: pair[0].z },
        end:   { x: pair[1].x, y: pair[1].y, z: pair[1].z },
        radius: .045, dashed: true, dashLength: .25, gapLength: .2,
        fromCap: false, toCap: false, color: "yellow"
      });
    });
    // In-scene labels for every nearby residue (not just polar partners) —
    // same devicePixelRatio supersampling as onStructureAtomClick's label,
    // so these read sharp too.
    var dpr = window.devicePixelRatio || 1;
    nearbyResnums.forEach(function (r) {
      var rep = nearbyReps[r];
      var lbl = structureViewer.addLabel((AA_3TO1[rep.resn] || rep.resn) + r, {
        position: { x: rep.x, y: rep.y, z: rep.z },
        backgroundColor: "#39424f", backgroundOpacity: .78,
        fontColor: "white", fontSize: 11 * dpr, padding: 3 * dpr, borderThickness: 0
      });
      if (lbl && lbl.sprite && dpr !== 1) lbl.sprite.scale.set(1 / dpr, 1 / dpr, 1);
    });
    structureViewer.render();

    var selfPill = neighborPillHtml(protein, resnum, targetAtoms[0].resn, true);
    var nearbyPills = nearbyResnums.length
      ? nearbyResnums.map(function (r) { return neighborPillHtml(protein, r, nearbyReps[r].resn, false); }).join("")
      : '<span class="sp-pill-none">none</span>';
    var polarPills = polarResnums.length
      ? polarResnums.map(function (r) { return neighborPillHtml(protein, r, polarResidues[r], false); }).join("")
      : '<span class="sp-pill-none">none found</span>';
    var text = document.getElementById("sp-neighbors");
    if (text) {
      text.innerHTML =
        "<div>" + selfPill + " is in proximity to (within " + NEIGHBOR_RADIUS + "Å): " + nearbyPills + "</div>" +
        "<div>Its side chain forms polar interactions with: " + polarPills + "</div>";
      text.style.display = "block";
    }
  }

  function closeStructurePanel() {
    var panel = document.getElementById("structure-panel");
    if (panel) panel.classList.remove("sp-open");
    document.body.classList.remove("structure-panel-open");
  }

  function setStructureStatus(msg) {
    var el = document.getElementById("sp-status");
    if (!el) return;
    el.textContent = msg || "";
    el.style.display = msg ? "block" : "none";
  }

  // Entry point: setupResidueClick() below calls this with whatever residue
  // the user just clicked on the Map.
  // Nearby Residues / All Variants / Compare all read atoms straight off
  // structureViewer's CURRENTLY LOADED model — clickable the instant the
  // panel opens, but a click before that model finishes loading (the very
  // first residue of a session, or a slow connection) would silently do
  // nothing (0 atoms match, 0 markers drawn, no error). Disabling them for
  // the loading window is simpler than making each one queue/retry.
  function setActionButtonsEnabled(enabled) {
    ["sp-neighbors-btn", "sp-allvariants-btn", "sp-compare-btn"].forEach(function (id) {
      var btn = document.getElementById(id);
      if (btn) btn.disabled = !enabled;
    });
  }

  function showResidueIn3D(protein, resnum, label) {
    buildStructurePanel();
    document.getElementById("structure-panel").classList.add("sp-open");
    document.body.classList.add("structure-panel-open");
    document.getElementById("sp-title").textContent = protein + " — " + (label || ("residue " + resnum));
    setStructureStatus("Loading " + protein + " structure…");
    setActionButtonsEnabled(false);

    load3Dmol()
      .then(function () { return fetchStructure(protein); })
      .then(function () {
        setStructureStatus(null);
        renderResidue(protein, resnum);
        setActionButtonsEnabled(true);
        updateCompareButton(protein, resnum);   // may re-disable Compare specifically (no analog at this residue)
      })
      .catch(function (err) {
        setStructureStatus("Couldn't load the 3D structure (" + err.message + ").");
      });
  }

  var NON_DOMAIN_KEYS = { linker: 1, nterm: 1, cterm: 1, none: 1 };
  function residueColor(protein, resnum) {
    if (structureColorMode === "domain") {
      var key = domainInfoFor(protein, resnum).key;
      return NON_DOMAIN_KEYS[key] ? STRUCTURE_NO_DOMAIN_GREY : domainColors()[key];
    }
    var score = conservationFor(protein, resnum);
    return (score != null && CONSERVATION_COLORS[score]) || STRUCTURE_NO_DOMAIN_GREY;
  }
  // Colors every atom's cartoon by the current mode, via a per-atom callback
  // (atom.resi is 3Dmol's own residue-number field off the parsed PDB).
  // setStyle REPLACES style for matched atoms rather than merging — calling
  // this on {} (all atoms) is what clears a previous highlight's stick/
  // sphere before highlightResidue() re-adds one, so the two always run
  // together (see setColorMode / highlightResidue's own callers). Also
  // (re-)establishes click-to-identify on every atom — setClickable has to
  // be re-applied whenever a new model is loaded (see renderResidue).
  function applyStructureStyle(protein) {
    structureViewer.setStyle({}, { cartoon: { colorfunc: function (atom) { return residueColor(protein, atom.resi); } } });
    structureViewer.setClickable({}, true, onStructureAtomClick);
  }
  // Layers a stick+sphere marker onto ONE residue via addStyle (adds a
  // representation without touching the cartoon applyStructureStyle() just
  // set — setStyle here would strip that residue's cartoon instead). "Jmol"
  // colorscheme for both (not a flat color) so heteroatoms in the side
  // chain read by element — N blue, O red, S yellow — the same convention
  // as PyMOL's default element coloring.
  function highlightResidue(resnum) {
    currentResnum = resnum;
    structureViewer.addStyle({ resi: resnum },
      { stick: { colorscheme: "Jmol", radius: .16 }, sphere: { colorscheme: "Jmol", scale: .22 } });
  }
  // Click-to-identify: any atom on the structure (not just the highlighted
  // residue) reports itself via a floating in-scene label, PyMOL-style.
  // Re-set on every render (see applyStructureStyle) since it targets the
  // CURRENTLY loaded model's atoms, not a fixed selection.
  function onStructureAtomClick(atom) {
    if (!structureViewer) return;
    var oneLetter = AA_3TO1[atom.resn] || atom.resn;
    structureViewer.removeAllLabels();
    // 3Dmol rasterizes a label's text onto a plain <canvas> at the raw
    // fontSize (CSS pixels), with no devicePixelRatio awareness — on a
    // retina display that texture then gets upscaled onto the (correctly
    // HiDPI-sized) WebGL canvas and looks soft, unlike the natively-drawn
    // cartoon/stick geometry. Supersample the texture (bigger fontSize/
    // padding) and scale the sprite back down by the same factor to land
    // back at the original apparent size, but sharper.
    var dpr = window.devicePixelRatio || 1;
    var label = structureViewer.addLabel(oneLetter + atom.resi, {
      position: { x: atom.x, y: atom.y, z: atom.z },
      backgroundColor: "#1c1f26", backgroundOpacity: .85,
      fontColor: "white", fontSize: 13 * dpr, padding: 4 * dpr, borderThickness: 0
    });
    if (label && label.sprite && dpr !== 1) label.sprite.scale.set(1 / dpr, 1 / dpr, 1);
    structureViewer.render();
  }

  function renderResidue(protein, resnum) {
    var el = document.getElementById("sp-viewer");
    if (!structureViewer) structureViewer = window.$3Dmol.createViewer(el, { backgroundColor: structureBg, antialias: true });
    if (loadedProtein !== protein) {
      structureViewer.removeAllModels();
      structureViewer.addModel(structureText[protein], "pdb");
      loadedProtein = protein;
    }
    structureViewer.removeAllLabels();   // clear any "click-to-identify" label from a prior residue
    structureViewer.removeAllShapes();   // clear any prior residue's "nearby" lines / "all variants" spheres
    showingNeighbors = false;
    showingAllVariants = false;
    var neighborsBtn = document.getElementById("sp-neighbors-btn");
    if (neighborsBtn) neighborsBtn.classList.remove("sp-fs-btn-active");
    var allVariantsBtn = document.getElementById("sp-allvariants-btn");
    if (allVariantsBtn) allVariantsBtn.classList.remove("sp-fs-btn-active");
    var neighborsText = document.getElementById("sp-neighbors");
    if (neighborsText) { neighborsText.innerHTML = ""; neighborsText.style.display = "none"; }
    applyStructureStyle(protein);
    highlightResidue(resnum);
    structureViewer.resize();   // panel may have just become visible; canvas size can be stale otherwise
    // No animation (duration 0 / omitted) — an earlier feature on this page
    // (the category-toggle redesign) hit real jank from animated transitions
    // depending on requestAnimationFrame timing; instant camera moves avoid
    // that whole class of bug here too.
    if (!hasFramedView) {
      // Very first click of the session — no established view to preserve,
      // so frame it fresh. zoomTo() alone frames just the clicked residue's
      // own atoms — tight enough that the surrounding fold (the actual
      // point of looking here) barely shows; pull back afterward for
      // context. Every later click, including a switch to the OTHER
      // protein, uses center() instead (below) so a zoom level the user set
      // up by scrolling survives — including across that switch.
      structureViewer.zoomTo({ resi: resnum });
      structureViewer.zoom(.4);
      hasFramedView = true;
    } else if (typeof structureViewer.center === "function") {
      structureViewer.center({ resi: resnum });
    } else {
      structureViewer.zoomTo({ resi: resnum });
    }
    structureViewer.render();

    var cap1 = document.getElementById("sp-pane-cap-1");
    if (cap1) cap1.textContent = protein + " — " + (AA_3TO1[targetAtoms0Resn(resnum)] || "") + resnum;
    updateCompareButton(protein, resnum);
    // Compare pane already open — keep it following the Map instead of
    // making the user re-click Compare for every new residue. If this
    // particular residue has no analogous position (an indel column),
    // just exit rather than show a stale/wrong comparison.
    if (compareMode) {
      var analog = analogousResidue(protein, resnum);
      if (analog != null) renderResidue2(otherProtein(protein), analog);
      else exitCompare();
    }
  }
  function targetAtoms0Resn(resnum) {
    var a = structureViewer.selectedAtoms({ resi: resnum })[0];
    return a ? a.resn : "";
  }
  function updateCompareButton(protein, resnum) {
    var btn = document.getElementById("sp-compare-btn");
    if (!btn) return;
    var other = otherProtein(protein);
    var analog = analogousResidue(protein, resnum);
    btn.textContent = "Compare to " + other;
    btn.disabled = analog == null;
    btn.title = analog == null ? "No analogous position in " + other + " for this residue" : "";
  }

  // Same click-to-identify behavior as onStructureAtomClick, but for
  // whichever viewer it's bound to — a factory instead of a second
  // hand-copied function, since pane 2 (Compare) needs its own instance
  // bound to structureViewer2, not the primary structureViewer.
  function makeAtomClickHandler(viewer) {
    return function (atom) {
      if (!viewer) return;
      var oneLetter = AA_3TO1[atom.resn] || atom.resn;
      viewer.removeAllLabels();
      var dpr = window.devicePixelRatio || 1;
      var label = viewer.addLabel(oneLetter + atom.resi, {
        position: { x: atom.x, y: atom.y, z: atom.z },
        backgroundColor: "#1c1f26", backgroundOpacity: .85,
        fontColor: "white", fontSize: 13 * dpr, padding: 4 * dpr, borderThickness: 0
      });
      if (label && label.sprite && dpr !== 1) label.sprite.scale.set(1 / dpr, 1 / dpr, 1);
      viewer.render();
    };
  }
  function applyStructureStyle2(protein) {
    structureViewer2.setStyle({}, { cartoon: { colorfunc: function (atom) { return residueColor(protein, atom.resi); } } });
    structureViewer2.setClickable({}, true, makeAtomClickHandler(structureViewer2));
  }
  function highlightResidue2(resnum) {
    currentResnum2 = resnum;
    structureViewer2.addStyle({ resi: resnum },
      { stick: { colorscheme: "Jmol", radius: .16 }, sphere: { colorscheme: "Jmol", scale: .22 } });
  }
  // Pane 2's camera is never independently framed — it just COPIES pane 1's
  // current view (setView, below). Both structures already share one
  // coordinate frame (the superposition align_structures.py computed), so
  // mirroring the raw camera state is what actually lands on "the same
  // viewing angle of the corresponding fold", more robust than re-deriving
  // a center point from the (sequence-based, not structure-based) residue
  // map — that map is only used to pick WHICH residue to highlight here.
  function renderResidue2(protein, resnum) {
    var el = document.getElementById("sp-viewer-2");
    if (!structureViewer2) structureViewer2 = window.$3Dmol.createViewer(el, { backgroundColor: structureBg, antialias: true });
    if (loadedProtein2 !== protein) {
      structureViewer2.removeAllModels();
      structureViewer2.addModel(structureText[protein], "pdb");
      loadedProtein2 = protein;
    }
    structureViewer2.removeAllLabels();
    applyStructureStyle2(protein);
    highlightResidue2(resnum);
    structureViewer2.resize();
    structureViewer2.setView(structureViewer.getView());
    structureViewer2.render();
    var resn = structureViewer2.selectedAtoms({ resi: resnum })[0];
    var cap2 = document.getElementById("sp-pane-cap-2");
    if (cap2) cap2.textContent = protein + " — " + (resn ? (AA_3TO1[resn.resn] || resn.resn) : "") + resnum;
  }

  function toggleCompare() {
    if (compareMode) { exitCompare(); return; }
    if (!loadedProtein || currentResnum == null) return;
    var other = otherProtein(loadedProtein);
    var otherResnum = analogousResidue(loadedProtein, currentResnum);
    if (otherResnum == null) return;   // button is disabled in this case (updateCompareButton) — belt and suspenders
    compareMode = true;
    document.getElementById("sp-compare-btn").classList.add("sp-fs-btn-active");
    document.getElementById("sp-viewer-pane-2").style.display = "block";
    resizeStructureViewer();   // pane 1 just shrank from 100% to 50% width
    setStructureStatus("Loading " + other + " structure…");
    load3Dmol()
      .then(function () { return fetchStructure(other); })
      .then(function () {
        setStructureStatus(null);
        renderResidue2(other, otherResnum);
      })
      .catch(function (err) {
        setStructureStatus("Couldn't load the comparison structure (" + err.message + ").");
        exitCompare();
      });
  }
  function exitCompare() {
    compareMode = false;
    var btn = document.getElementById("sp-compare-btn");
    if (btn) btn.classList.remove("sp-fs-btn-active");
    var pane2 = document.getElementById("sp-viewer-pane-2");
    if (pane2) pane2.style.display = "none";
    resizeStructureViewer();   // pane 1 back to the full row width
  }

  // Both viewers already share one coordinate frame (superposition), so
  // "synced rotation" is just: on any mouse/touch/wheel interaction with
  // EITHER viewer's canvas, copy its current view onto the other. 3Dmol
  // updates the camera synchronously inside ITS OWN listener on the same
  // element (attached first, during createViewer) — by the time this
  // listener runs, getView() already reflects that interaction.
  // Deliberately NOT requestAnimationFrame-throttled: this page already
  // learned that lesson once (the category-toggle redesign hit real jank
  // from rAF-dependent timing) — a plain synchronous copy on every event
  // costs one extra render() per input event, no worse than what 3Dmol's
  // OWN drag handler is already doing on that same element, and it can't
  // go stale or get silently dropped the way an rAF callback can. No
  // reentrancy guard needed either: copying a view onto the OTHER canvas
  // via setView()+render() doesn't itself dispatch a mouse/wheel DOM event,
  // so there's no ping-pong loop to guard against.
  function syncFromPane1() {
    if (!compareMode || !structureViewer || !structureViewer2) return;
    structureViewer2.setView(structureViewer.getView());
    structureViewer2.render();
  }
  function syncFromPane2() {
    if (!compareMode || !structureViewer || !structureViewer2) return;
    structureViewer.setView(structureViewer2.getView());
    structureViewer.render();
  }
  function wireCompareSync() {
    var v1 = document.getElementById("sp-viewer"), v2 = document.getElementById("sp-viewer-2");
    ["mousemove", "wheel", "touchmove"].forEach(function (evt) {
      v1.addEventListener(evt, syncFromPane1, { passive: true });
      v2.addEventListener(evt, syncFromPane2, { passive: true });
    });
  }

  // Click ANY residue box on the Map (curated variant or not — this is
  // separate from the .vtick/.vlabel popup click handler in setupPopup(),
  // which lives on different overlay elements entirely, so the two never
  // fire for the same click) to open/update the 3D panel.
  function setupResidueClick(table) {
    if (table.__resClick) return;
    table.__resClick = true;
    table.addEventListener("click", function (e) {
      var cell = e.target.closest ? e.target.closest("td[data-pos]") : null;
      if (!cell) return;
      var protein = cell.dataset.protein;
      var resnum = parseInt(cell.dataset.resnum, 10);
      if (!protein || !resnum) return;
      showResidueIn3D(protein, resnum, cell.dataset.pos);
    });
  }

  // toggleBox is the ONE actual toggle-list element — never duplicated.
  // On Map it's a fixed top-right box; on Table it's re-parented into a
  // permanently-visible row under the other filter controls (see
  // dockCategoryToggles/setView below). Same checkboxes, same listener,
  // wherever it currently lives.
  var toggleBox = null;
  function buildToggles() {
    if (document.getElementById("variant-toggles")) return;
    var box = document.createElement("div");
    box.id = "variant-toggles";
    var html = '<div class="vt-title">Show variants</div>';
    Object.keys(CFG).forEach(function (cat) {
      var on = CFG[cat].on !== false;
      if (!on) document.body.classList.add("hide-cat-" + cat);   // default state
      html += '<label><input type="checkbox" ' + (on ? "checked" : "") +
              ' data-cat="' + cat + '">' +
              '<span class="vt-name" style="--vt-c:' + CFG[cat].color + '">' +
              CFG[cat].legend + '</span></label>';
    });
    box.innerHTML = html;
    box.addEventListener("change", function (e) {
      var cb = e.target;
      if (cb.tagName !== "INPUT") return;
      document.body.classList.toggle("hide-cat-" + cb.getAttribute("data-cat"), !cb.checked);
      // same category toggles also gate table-view rows (they share the
      // data-category attribute + hide-cat-* CSS rules); just the visible
      // count needs an explicit refresh since it's plain text, not CSS.
      var tv = document.getElementById("table-view");
      if (tv && tv.__wired) updateTableCount();
    });
    document.body.appendChild(box);
    toggleBox = box;
  }
  // Re-parents the ONE toggleBox between its two homes. Table view builds a
  // fresh #vtbl-cat-row on every rebuild (sort clicks, gnomAD injection), so
  // this has to be callable repeatedly, not just once on the first switch.
  // Always fully visible in both homes — no button/dropdown, no hidden state.
  function dockCategoryToggles(view) {
    if (!toggleBox) return;
    if (view === "table") {
      var host = document.getElementById("table-view");
      var row = host && host.querySelector("#vtbl-cat-row");
      if (!row) return;
      row.appendChild(toggleBox);
      toggleBox.classList.add("vt-inline");
    } else {
      if (toggleBox.parentElement !== document.body) document.body.appendChild(toggleBox);
      toggleBox.classList.remove("vt-inline");
    }
  }
  // Adds ONE checkbox to the already-built toggle box, for a category that
  // only becomes known after the initial render (gnomADmis — discovered
  // once the live gnomAD fetch resolves). The box's existing "change"
  // listener is delegated (checks e.target), so this new row is covered by
  // it automatically — no separate listener needed, and any checkbox state
  // the user has already set on OTHER categories is left untouched.
  function addToggleFor(cat) {
    var box = document.getElementById("variant-toggles");
    if (!box || box.querySelector('input[data-cat="' + cat + '"]')) return;
    var on = CFG[cat].on !== false;
    if (!on) document.body.classList.add("hide-cat-" + cat);
    var label = document.createElement("label");
    label.innerHTML = '<input type="checkbox" ' + (on ? "checked" : "") + ' data-cat="' + cat + '">' +
                       '<span class="vt-name" style="--vt-c:' + CFG[cat].color + '">' + CFG[cat].legend + '</span>';
    box.appendChild(label);
  }

  // ---- Map / Table view switcher (top-left, mirrors the toggle box) ----
  // "Table" swaps the big aligned-sequence map for a sortable/filterable list
  // of every variant. The category toggles (top-right) are NOT duplicated —
  // table rows carry the same data-category attribute as ticks/labels, so the
  // existing body.hide-cat-* CSS rules already show/hide them for free.
  function buildViewSwitcher() {
    if (document.getElementById("view-switcher")) return;
    var box = document.createElement("div");
    box.id = "view-switcher";
    box.innerHTML =
      '<button type="button" class="vs-btn vs-active" data-view="map">Map</button>' +
      '<button type="button" class="vs-btn" data-view="table">Table</button>';
    box.addEventListener("click", function (e) {
      var btn = e.target.closest && e.target.closest("button[data-view]");
      if (btn) setView(btn.getAttribute("data-view"));
    });
    document.body.appendChild(box);
  }
  var currentView = "map";
  function setView(view) {
    if (view === currentView) return;
    var isTable = view === "table";
    var mapTable = document.getElementById("samd9-table");
    var tv = document.getElementById("table-view");
    var strip = document.getElementById("bottom-legend-strip");
    var bar = document.getElementById("bottom-legend-bar");
    if (mapTable) mapTable.style.display = isTable ? "none" : "";
    if (strip) strip.style.display = isTable ? "none" : "";
    // bar's own static HTML sets "display: flex" inline (it lays out its three
    // legend-blocks side by side) — restoring "" here would drop back to a
    // <div>'s default block display instead, stacking the conservation scale
    // above the other legend blocks instead of beside them. Must name "flex"
    // explicitly, not just clear the property.
    if (bar) bar.style.display = isTable ? "none" : "flex";
    if (tv) tv.style.display = isTable ? "block" : "none";
    document.querySelectorAll("#view-switcher .vs-btn").forEach(function (b) {
      b.classList.toggle("vs-active", b.getAttribute("data-view") === view);
    });
    dockCategoryToggles(view);
    currentView = view;
    // run() no-ops while hidden (see its own comment), so geometry can go
    // stale while Table was showing (e.g. the window was resized). Force one
    // fresh, correct pass now that the Map is visible again. Deferred one
    // tick (not immediate) so the display:"" swap above gets to paint before
    // run()'s own heavy reflow (it walks every tick/label + domain-outline
    // cell) runs on the main thread; run() itself already no-ops if the user
    // has switched away again by the time this fires, so a rapid
    // Map->Table->Map isn't at risk of doing wasted work.
    if (!isTable) setTimeout(run, 0);
  }

  // ---- Table view: every variant as a sortable, filterable list ----
  // `width` becomes a <colgroup><col> hint (buildTableView, below), pinning
  // every column to a fixed size so a long value (e.g. "SIR2–P-loop NTPase
  // linker") can never overflow into its neighbor the way table-layout:auto
  // let it — auto was dividing leftover space EQUALLY across columns
  // regardless of content, so Domain and Conservation ended up the same
  // width despite needing very different amounts of room. The last column
  // (Source) has no width — it absorbs whatever's left.
  var TABLE_COLUMNS = [
    { key: "protein",      label: "Protein",             width: 95 },
    { key: "variant",      label: "Variant",             width: 130 },
    { key: "conservation", label: 'Conservation<br><span class="vtbl-th-sub">score</span>', width: 100 },
    { key: "domain",       label: "Domain",              width: 260 },
    { key: "site",         label: 'Nearest functional<br><span class="vtbl-th-sub">site</span>', width: 170 },
    { key: "category",     label: "Category",           width: 210 },
    { key: "method",       label: "Method",              width: 150 },
    { key: "gnomad",       label: 'gnomAD<br><span class="vtbl-th-sub">frequency</span>',   width: 120 },
    { key: "hom",          label: 'gnomAD<br><span class="vtbl-th-sub">homozygotes</span>', width: 130 },
    { key: "source",       label: "Source" }
  ];
  var EXPORT_HEADERS = ["Protein", "Variant", "Conservation", "Domain", "Nearest functional site", "Distance (Å)", "Category", "Method", "gnomAD frequency", "gnomAD homozygotes", "Source", "PMID"];
  var SORT_VAL = {
    protein:      function (r) { return r.protein; },
    variant:      function (r) { return r.label.toLowerCase(); },
    conservation: function (r) { return r.conservation == null ? -1 : r.conservation; },
    domain:       function (r) { return DOMAIN_FILTER_ORDER.indexOf(r.domain); },
    // No site data (residue not found in the AlphaFold model — an
    // alignment-gap position) sorts to the end in both directions: always
    // Infinity, never negative, so ascending doesn't put "unknown" first.
    site:         function (r) { return r.site ? r.site.distance : Infinity; },
    category:     function (r) { return (CFG[r.category] || CFG.Other).legend.toLowerCase(); },
    method:       function (r) { return (r.method || "").toLowerCase(); },
    gnomad:       function (r) { return r.gnomadAF == null ? -1 : r.gnomadAF; },
    hom:          function (r) { return r.gnomadHom == null ? -1 : r.gnomadHom; },
    source:       function (r) { return (r.sourceTitle || "").toLowerCase(); }
  };
  var tableSort = { key: "protein", dir: 1 };
  // Both Domain and Conservation are multi-select (checkbox panels, not a
  // single dropdown) — same mechanics for both, so MULTI_FILTERS/tableFilter
  // .multi drive both through one set of shared functions instead of two
  // near-duplicate ones. Each starts holding every one of its own options
  // (nothing filtered out) rather than an "all" sentinel.
  var CONSERVATION_SCORES = ["9", "8", "7", "6", "5", "4", "3", "2", "1"];
  var MULTI_FILTERS = {
    domain: { options: DOMAIN_FILTER_ORDER, label: "Domain" },
    conservation: { options: CONSERVATION_SCORES, label: "Conservation" }
  };
  // Search text + protein pick are table-view-only state too, independent of
  // the shared category toggles — search is a plain substring match against
  // each row's own rendered text (already covers protein/variant/domain/
  // category/method/gnomAD/source in one go, no separate index).
  var tableFilter = {
    search: "", protein: "all",
    multi: { domain: new Set(DOMAIN_FILTER_ORDER), conservation: new Set(CONSERVATION_SCORES) }
  };

  // Same present/absent/loading logic as gnomadRowFor (the popup), but split
  // into a display string + numeric allele frequency + homozygote count so
  // those columns can be sorted, not just read.
  function gnomadCellFor(protein, label, staticGnomad) {
    if (GNOMAD_LIVE.ready && isSimpleLabel(label)) {
      var hit = GNOMAD_LIVE.byLabel[protein] && GNOMAD_LIVE.byLabel[protein][label];
      // gnomAD's API returns homozygote_count as null for some present
      // variants (not computed/redacted, seemingly independent of how rare
      // the variant is) — treat that as 0, not "unknown": every variant we
      // can confirm is PRESENT should show a real number here, so "—" means
      // only "not present / status unknown", never "present but who knows".
      if (hit && hit.ac > 0) return { text: hit.af.toExponential(2), af: hit.af, hom: hit.hom == null ? 0 : hit.hom };
      return { text: "Not present", af: 0, hom: null };
    }
    if (staticGnomad != null) {
      var g = (typeof staticGnomad === "object")
        ? (staticGnomad.present ? "Yes" + (staticGnomad.maf ? " (MAF " + staticGnomad.maf + ")" : "") : "No")
        : staticGnomad;
      return { text: g, af: null, hom: null };
    }
    if (GNOMAD_LIVE.failed) return { text: "—", af: null, hom: null };
    return { text: isSimpleLabel(label) ? "Loading…" : "—", af: null, hom: null };
  }

  // Rebuilds the whole table from the current DATA snapshot. Called after
  // init (DATA loaded) and after gnomAD-missense injection (DATA grew) — NOT
  // on every resize/run(), since rebuilding thousands of rows on resize would
  // be wasteful and gains nothing (the table has no geometry to recompute).
  // `lastRows` (module-level) is the export buttons' source: they read
  // whichever of these rows correspond to currently-visible <tr>s.
  var lastRows = [];
  function buildTableView() {
    var host = document.getElementById("table-view");
    if (!host) return;

    var rows = DATA.map(function (v) {
      var d = DETAILS[v.protein + ":" + v.label] || DETAILS[v.label];
      var paper = (d && d.paper && PAPERS[d.paper]) || null;
      var gc = gnomadCellFor(v.protein, v.label, d && d.gnomad);
      var dom = domainInfoFor(v.protein, v.residue);
      return {
        protein: v.protein, residue: v.residue, label: v.label, category: v.category,
        domain: dom.name, domainKey: dom.key,
        conservation: conservationFor(v.protein, v.residue),
        site: siteFor(v.protein, v.residue),
        method: d && d.method,
        gnomadAF: gc.af, gnomadText: gc.text, gnomadHom: gc.hom,
        sourceTitle: paper && paper.title, sourceUrl: paper && paper.url, pmid: paper && paper.pmid
      };
    });

    var valFn = SORT_VAL[tableSort.key] || SORT_VAL.protein;
    rows.sort(function (a, b) {
      var fa = valFn(a), fb = valFn(b);
      var cmp = fa < fb ? -1 : fa > fb ? 1 : 0;
      if (!cmp) cmp = (a.protein < b.protein ? -1 : a.protein > b.protein ? 1 : 0) || (a.residue - b.residue);
      return cmp * tableSort.dir;
    });
    lastRows = rows;

    var colgroupHtml = "<colgroup>" + TABLE_COLUMNS.map(function (c) {
      return c.width ? '<col style="width:' + c.width + 'px">' : "<col>";
    }).join("") + "</colgroup>";

    var theadHtml = "<tr>" + TABLE_COLUMNS.map(function (c) {
      var arrow = tableSort.key === c.key ? (tableSort.dir === 1 ? " ▲" : " ▼") : "";
      // c.label is a hardcoded constant (never variant data), so the "<br>"
      // a couple of labels carry (to wrap "gnomAD" onto its own line — those
      // two headers were overflowing even on a large screen) can go through
      // unescaped here.
      return '<th data-sort="' + c.key + '">' + c.label + arrow + "</th>";
    }).join("") + "</tr>";

    var dcolors = domainColors();
    var bodyHtml = rows.map(function (r) {
      var color = (CFG[r.category] || CFG.Other).color;
      // The toggle box keeps "(unannotated)" — useful context there, next to
      // GoF/LoF/Somatic/etc in a legend the reader sees once. Repeated down
      // a table column hundreds of times it's just noise, so the Category
      // badge drops it; nothing else about the category changes.
      var legend = (CFG[r.category] || CFG.Other).legend.replace(" (unannotated)", "");
      var src = "—";
      if (r.sourceTitle) {
        src = (r.sourceUrl && r.sourceUrl !== "pending")
          ? '<a href="' + esc(r.sourceUrl) + '" target="_blank" rel="noopener" title="' + esc(r.sourceTitle) + '">' + esc(r.sourceTitle) + "</a>"
          : esc(r.sourceTitle);
        if (r.pmid && r.pmid !== "pending") src += '<div class="vtbl-pmid">PMID: ' + esc(r.pmid) + "</div>";
      }
      // Zero and "unknown" were both landing on screen as visually distinct
      // (0 vs —), which read as an inconsistency rather than two different
      // real states. Simplify: only a CONFIRMED positive homozygote count
      // gets a number; zero and unknown both just show "—" on screen. (The
      // export still writes the real 0 — that distinction is worth keeping
      // in a data file, just not worth the on-screen confusion.)
      var homCell = r.gnomadHom > 0
        ? '<span class="vtbl-hom-pos">' + r.gnomadHom + "</span>"
        : "—";
      var consCell = r.conservation != null
        ? '<span class="cs-box vtbl-cons Score' + r.conservation + '">' + r.conservation + "</span>"
        : "—";
      var siteCell = r.site
        ? r.site.distance + " Å<div class=\"vtbl-pmid\">" + esc(r.site.site) + "</div>"
        : "—";
      return '<tr data-category="' + esc(r.category) + '" data-protein="' + esc(r.protein) + '" data-domain="' + esc(r.domain) +
        '" data-conservation="' + (r.conservation == null ? "" : r.conservation) + '">' +
        "<td>" + esc(r.protein) + "</td>" +
        '<td class="vtbl-mono">' + esc(r.label) + "</td>" +
        "<td>" + consCell + "</td>" +
        '<td><span class="vtbl-dom" style="--dm-c:' + dcolors[r.domainKey] + '">' + esc(r.domain) + "</span></td>" +
        "<td>" + siteCell + "</td>" +
        '<td class="vtbl-left"><span class="vtbl-cat" style="--vt-c:' + color + '">' + esc(legend) + "</span></td>" +
        // <wbr> (a soft break hint) let the browser's greedy line-fill still
        // choose to break at the LATER space instead — "Viral infection/flow"
        // / "cytometry" — since that packed more characters onto line 1. A
        // hard <br> forces the split to always land right after the "/".
        "<td>" + (r.method ? esc(r.method).replace(/\//g, "/<br>") : "—") + "</td>" +
        "<td>" + esc(r.gnomadText) + "</td>" +
        "<td>" + homCell + "</td>" +
        '<td class="vtbl-source">' + src + "</td>" +
        "</tr>";
    }).join("");

    var proteins = ["all", "SAMD9", "SAMD9L"];
    var proteinHtml = proteins.map(function (p) {
      var active = tableFilter.protein === p ? " vs-active" : "";
      return '<button type="button" class="vs-btn vtbl-pbtn' + active + '" data-protein="' + p + '">' +
             (p === "all" ? "All" : p) + "</button>";
    }).join("");
    var domainPanelHtml = multiSelectHtml("domain");
    var conservationPanelHtml = multiSelectHtml("conservation");

    host.innerHTML =
      '<div class="vtbl-head">' +
        '<div class="vtbl-controls">' +
          '<input type="search" id="vtbl-search" placeholder="Search variants…" value="' + esc(tableFilter.search) + '">' +
          '<div class="vtbl-protein-filter">' + proteinHtml + "</div>" +
          domainPanelHtml +
          conservationPanelHtml +
          '<button type="button" class="vtbl-export" id="vtbl-export-csv">Export CSV</button>' +
          '<button type="button" class="vtbl-export" id="vtbl-export-xlsx">Export Excel</button>' +
        "</div>" +
        '<span id="vtbl-count"></span>' +
      "</div>" +
      '<div class="vtbl-cat-row" id="vtbl-cat-row"></div>' +
      '<div class="vtbl-scroll"><table class="vtbl">' + colgroupHtml + "<thead>" + theadHtml + "</thead><tbody>" + bodyHtml + "</tbody></table></div>";
    applyTableFilters();   // re-apply protein/domain(s)/conservation(s)/search to the freshly-built rows
    // The whole host got fully replaced above, which orphans the category
    // checkboxes if they were docked in the old #vtbl-cat-row (sort clicks,
    // gnomAD injection both call buildTableView() again) — re-home them in
    // the fresh row. Only when Table is actually the active view: this also
    // runs from init() before the user has ever switched views, and toggleBox
    // must stay in its Map-mode fixed position until they actually do.
    if (currentView === "table") dockCategoryToggles("table");

    if (!host.__wired) {
      host.__wired = true;
      host.addEventListener("click", function (e) {
        var th = e.target.closest && e.target.closest("th[data-sort]");
        if (th) {
          var key = th.getAttribute("data-sort");
          if (tableSort.key === key) tableSort.dir *= -1; else { tableSort.key = key; tableSort.dir = 1; }
          buildTableView();
          return;
        }
        var pbtn = e.target.closest && e.target.closest(".vtbl-pbtn");
        if (pbtn) {
          tableFilter.protein = pbtn.getAttribute("data-protein");
          host.querySelectorAll(".vtbl-pbtn").forEach(function (b) {
            b.classList.toggle("vs-active", b === pbtn);
          });
          applyTableFilters();
          return;
        }
        var mbtn = e.target.closest && e.target.closest(".vtbl-mbtn");
        if (mbtn) {
          var mt = mbtn.getAttribute("data-mtarget");
          var panel = document.getElementById("vtbl-" + mt + "-panel");
          if (panel) panel.toggleAttribute("hidden");
          return;
        }
        var daction = e.target.closest && e.target.closest("[data-daction]");
        if (daction) {
          var mkey = daction.getAttribute("data-mtarget");
          var toAll = daction.getAttribute("data-daction") === "all";
          tableFilter.multi[mkey] = new Set(toAll ? MULTI_FILTERS[mkey].options : []);
          host.querySelectorAll('.vtbl-mcb[data-mtarget="' + mkey + '"]').forEach(function (cb) { cb.checked = toAll; });
          document.getElementById("vtbl-" + mkey + "-btn").textContent = multiSelectLabel(mkey);
          applyTableFilters();
          return;
        }
        if (e.target.id === "vtbl-export-csv") { exportCSV(); return; }
        if (e.target.id === "vtbl-export-xlsx") { exportExcel(); return; }
        // Clicking anywhere else while a filter panel is open closes it — but
        // not a click ON that panel/its checkboxes/its own button (each of
        // those is handled above, or is a checkbox toggle via the "change"
        // listener below, and none of those paths should also close it).
        Object.keys(MULTI_FILTERS).forEach(function (key) {
          var panel = document.getElementById("vtbl-" + key + "-panel");
          if (panel && !panel.hasAttribute("hidden") && !e.target.closest('[data-mtarget="' + key + '"]')) {
            panel.setAttribute("hidden", "");
          }
        });
      });
      host.addEventListener("change", function (e) {
        if (e.target.classList && e.target.classList.contains("vtbl-mcb")) {
          var mkey2 = e.target.getAttribute("data-mtarget");
          if (e.target.checked) tableFilter.multi[mkey2].add(e.target.value);
          else tableFilter.multi[mkey2].delete(e.target.value);
          document.getElementById("vtbl-" + mkey2 + "-btn").textContent = multiSelectLabel(mkey2);
          applyTableFilters();
        }
      });
      var searchTimer;
      host.addEventListener("input", function (e) {
        if (e.target.id !== "vtbl-search") return;
        var val = e.target.value;
        clearTimeout(searchTimer);
        searchTimer = setTimeout(function () {
          tableFilter.search = val.trim().toLowerCase();
          applyTableFilters();
        }, 150);
      });
    }
  }
  // Builds one multi-select filter (a button that opens a checkbox panel) —
  // shared by Domain and Conservation, the two filters where more than one
  // value can be picked at once (unlike the plain protein buttons).
  function multiSelectHtml(key) {
    var cfg = MULTI_FILTERS[key];
    var set = tableFilter.multi[key];
    var cbHtml = cfg.options.map(function (name) {
      var checked = set.has(name) ? " checked" : "";
      return '<label><input type="checkbox" class="vtbl-mcb" data-mtarget="' + key + '" value="' + esc(name) + '"' + checked + '><span>' + esc(name) + "</span></label>";
    }).join("");
    return '<div class="vtbl-domain-wrap" data-mtarget="' + key + '">' +
      '<button type="button" class="vtbl-export vtbl-mbtn" data-mtarget="' + key + '" id="vtbl-' + key + '-btn">' + esc(multiSelectLabel(key)) + "</button>" +
      '<div id="vtbl-' + key + '-panel" class="vtbl-domain-panel" hidden>' +
        '<div class="vtbl-domain-actions">' +
          '<button type="button" data-mtarget="' + key + '" data-daction="all">Select all</button>' +
          '<button type="button" data-mtarget="' + key + '" data-daction="none">Clear</button>' +
        "</div>" + cbHtml +
      "</div>" +
    "</div>";
  }
  function multiSelectLabel(key) {
    var cfg = MULTI_FILTERS[key];
    var n = tableFilter.multi[key].size, total = cfg.options.length;
    if (n === total) return cfg.label + ": All ▾";
    if (n === 0) return cfg.label + ": None ▾";
    return cfg.label + ": " + n + " selected ▾";
  }
  // Applies the protein/domain(s)/conservation(s) pick + search text on top
  // of whatever the category toggles already did via CSS (their !important
  // rule always wins, so we never fight it — a category-hidden row just
  // stays hidden regardless of what we set here). Cheap enough to run on
  // every keystroke: substring match against the row's own rendered text,
  // no separate search index.
  function applyTableFilters() {
    var host = document.getElementById("table-view");
    if (!host) return;
    var rows = host.querySelectorAll("tbody tr");
    rows.forEach(function (tr) {
      var proteinOk = tableFilter.protein === "all" || tr.getAttribute("data-protein") === tableFilter.protein;
      var domainOk = tableFilter.multi.domain.has(tr.getAttribute("data-domain"));
      var consAttr = tr.getAttribute("data-conservation");
      var consOk = consAttr === "" || tableFilter.multi.conservation.has(consAttr);
      var searchOk = !tableFilter.search || tr.textContent.toLowerCase().indexOf(tableFilter.search) !== -1;
      tr.style.display = (proteinOk && domainOk && consOk && searchOk) ? "" : "none";
    });
    updateTableCount();
  }
  function updateTableCount() {
    var el = document.getElementById("vtbl-count");
    if (!el) return;
    var rows = document.querySelectorAll("#table-view tbody tr");
    var total = rows.length, shown = 0;
    rows.forEach(function (tr) { if (getComputedStyle(tr).display !== "none") shown++; });
    el.textContent = shown.toLocaleString() + " of " + total.toLocaleString() + " variants shown";
  }

  // ---- Table export (CSV / Excel) — exports exactly what's currently
  // visible (category toggles + protein/domain filter + search), in the
  // current sort order, since that's "what you're looking at" in the table.
  function exportRowValues(r) {
    return [
      r.protein, r.label, r.conservation != null ? r.conservation : "", r.domain,
      r.site ? r.site.site : "", r.site ? r.site.distance : "",
      (CFG[r.category] || CFG.Other).legend.replace(" (unannotated)", ""),
      r.method || "", r.gnomadText, r.gnomadHom != null ? r.gnomadHom : "",
      r.sourceTitle || "", r.pmid || ""
    ];
  }
  function visibleExportRows() {
    var host = document.getElementById("table-view");
    if (!host) return [];
    var trs = host.querySelectorAll("tbody tr");
    var out = [];
    trs.forEach(function (tr, i) {
      if (getComputedStyle(tr).display !== "none" && lastRows[i]) out.push(lastRows[i]);
    });
    return out;
  }
  function downloadBlob(content, filename, mime) {
    var blob = new Blob([content], { type: mime });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }
  function csvEscape(v) {
    var s = String(v == null ? "" : v);
    return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function exportCSV() {
    var lines = [EXPORT_HEADERS.map(csvEscape).join(",")];
    visibleExportRows().forEach(function (r) { lines.push(exportRowValues(r).map(csvEscape).join(",")); });
    downloadBlob(lines.join("\r\n"), "samd9_samd9l_variants.csv", "text/csv;charset=utf-8;");
  }
  // A genuine Excel file with no external library: SpreadsheetML (Excel's
  // 2003 XML format) is a plain, documented XML schema Excel opens natively
  // and warning-free — unlike the common "save an HTML table as .xls" trick,
  // which triggers a format-mismatch warning. Colors category/domain cells to
  // match the on-page badges (lighten() reused from the domain-outline code
  // above) and highlights homozygote-positive counts, same as the table.
  function toExcelXML(rows) {
    var dcolors = domainColors();
    function cell(value, styleId, type) {
      var attr = styleId ? ' ss:StyleID="' + styleId + '"' : "";
      return "<Cell" + attr + '><Data ss:Type="' + (type || "String") + '">' + esc(value) + "</Data></Cell>";
    }
    var styles = [
      '<Style ss:ID="sHeader"><Font ss:Bold="1"/><Interior ss:Color="#F0F0F0" ss:Pattern="Solid"/></Style>',
      '<Style ss:ID="sHomPos"><Font ss:Color="#1E9E5B" ss:Bold="1"/></Style>'
    ];
    Object.keys(CFG).forEach(function (cat) {
      styles.push('<Style ss:ID="cat_' + cat + '"><Interior ss:Color="' + lighten(CFG[cat].color, 0.84) + '" ss:Pattern="Solid"/></Style>');
    });
    Object.keys(dcolors).forEach(function (key) {
      styles.push('<Style ss:ID="dom_' + key + '"><Interior ss:Color="' + lighten(dcolors[key], 0.82) + '" ss:Pattern="Solid"/></Style>');
    });
    var headerRow = "<Row>" + EXPORT_HEADERS.map(function (h) { return cell(h, "sHeader"); }).join("") + "</Row>";
    var bodyRows = rows.map(function (r) {
      var v = exportRowValues(r);
      return "<Row>" +
        cell(v[0]) + cell(v[1]) +
        cell(v[2], null, r.conservation != null ? "Number" : "String") +
        cell(v[3], "dom_" + r.domainKey) +
        cell(v[4]) + cell(v[5], null, r.site ? "Number" : "String") +
        cell(v[6], "cat_" + r.category) +
        cell(v[7]) + cell(v[8]) +
        cell(v[9], r.gnomadHom > 0 ? "sHomPos" : null, r.gnomadHom != null ? "Number" : "String") +
        cell(v[10]) + cell(v[11]) +
        "</Row>";
    }).join("");
    return '<?xml version="1.0"?><?mso-application progid="Excel.Sheet"?>' +
      '<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet" ' +
      'xmlns:o="urn:schemas-microsoft-com:office:office" ' +
      'xmlns:x="urn:schemas-microsoft-com:office:excel" ' +
      'xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet">' +
      "<Styles>" + styles.join("") + "</Styles>" +
      '<Worksheet ss:Name="Variants"><Table>' + headerRow + bodyRows + "</Table></Worksheet>" +
      "</Workbook>";
  }
  function exportExcel() {
    downloadBlob(toExcelXML(visibleExportRows()), "samd9_samd9l_variants.xls", "application/vnd.ms-excel");
  }

  // "Last updated" — fetched live from GitHub's own commit history, purely
  // client-side. No CI, no deploy step, nothing that can conflict with GitHub
  // Pages' own deployment — if the fetch fails (offline, rate-limited, repo
  // renamed) the static fallback date already in the HTML is left as-is.
  function updateLastUpdatedDate() {
    var el = document.querySelector("#last-updated .lu-date");
    if (!el) return;
    fetch("https://api.github.com/repos/pmitev93/SAMD9_L_Map/commits/main")
      .then(function (r) { return r.ok ? r.json() : Promise.reject(); })
      .then(function (d) {
        var iso = d.commit && d.commit.committer && d.commit.committer.date;
        if (!iso) return;
        var dt = new Date(iso);
        var months = ["January", "February", "March", "April", "May", "June",
                       "July", "August", "September", "October", "November", "December"];
        el.textContent = dt.getUTCDate() + " " + months[dt.getUTCMonth()] + " " + dt.getUTCFullYear();
      })
      .catch(function () { /* keep the static fallback already in the HTML */ });
  }

  var t;
  function schedule() { clearTimeout(t); t = setTimeout(run, 120); }

  // Load the editable data files fresh each time so "edit a data_*.js + refresh"
  // always shows your changes (browsers otherwise cache <script src>). On a
  // double-clicked file:// page we skip the query (not needed, and keeps it simple).
  function boot() {
    updateLastUpdatedDate();
    loadGnomadLive();   // fires in the background; popups just check GNOMAD_LIVE whenever opened
    var files = ["data_variants.js", "data_overrides.js", "data_papers.js", "data_residue_map.js", "data_site_distances.js"];
    var bust = location.protocol === "file:" ? "" : ("?t=" + Date.now());
    var left = files.length;
    files.forEach(function (f) {
      var s = document.createElement("script");
      s.src = f + bust;
      s.onload = s.onerror = function () { if (--left === 0) init(); };
      document.head.appendChild(s);
    });
  }
  if (document.readyState === "loading")
    document.addEventListener("DOMContentLoaded", boot);
  else boot();
  window.addEventListener("resize", schedule);
})();
