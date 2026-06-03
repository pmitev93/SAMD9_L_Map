/* Data-driven variant renderer for the SAMD9/SAMD9L map.
   Reads embedded variant-data (+ optional variant-overrides), removes the old
   hand-placed markers, and draws colored ticks + lane-stacked labels above
   SAMD9L / below SAMD9. Also builds the show/hide category toggles. */
(function () {
  "use strict";

  // ---- Toggle/legend text: edit the "legend" strings below to rename a category.
  var CFG = {
    GoF:    { color: "#FF0000", label: true,  on: true,  legend: "Gain-of-function" },
    LoF:    { color: "#1f56bc", label: true,  on: true,  legend: "Loss-of-function" },
    gnomAD: { color: "#73d73c", label: false, on: false, legend: "gnomAD (truncating)" },
    Other:  { color: "#000000", label: true,  on: true,  legend: "Other (somatic / NoF)" }
  };
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

  // ---- data: read from the data_*.js globals (loaded via <script src>, which
  //      works even on a double-clicked file). Edit a data_*.js + refresh = live.
  var DATA = [], OV = {}, PAPERS = {}, DETAILS = {};
  function ovFor(v) { return OV[v.label] || OV[v.protein + ":" + v.label] || null; }

  function init() {
    DATA    = window.VARIANTS || [];
    OV      = window.VARIANT_OVERRIDES || {};
    PAPERS  = window.PAPERS || {};
    DETAILS = window.VARIANT_DETAILS || {};
    run();
    setupPopup();
  }

  function run() {
    var table = document.querySelector("table");
    if (!table) return;

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

    // 4. residue -> cell map (cumulative per protein, document order)
    var maps = { SAMD9: {}, SAMD9L: {} };
    var cellRow = new Map();
    var counters = { SAMD9: 0, SAMD9L: 0 };
    seqRows.forEach(function (sr) {
      var cells = sr.row.children;
      for (var i = 1; i < cells.length; i++) {
        var txt = (cells[i].textContent || "").trim();
        if (/^[A-Z]$/.test(txt)) {
          counters[sr.protein]++;
          maps[sr.protein][counters[sr.protein]] = cells[i];
          cellRow.set(cells[i], sr);
          cells[i].dataset.pos = txt + counters[sr.protein];   // e.g. "K133"
          cells[i].dataset.side = sr.side;                     // top / bottom
        }
      }
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

      var lanes = [], maxLane = -1;
      list.forEach(function (v) {
        var cfg = CFG[v.category] || CFG.Other;
        v._color = cfg.color;
        if (cfg.label) {
          var w = estW(v.label), leftEdge = v._x - w / 2, l = 0;
          while (l < lanes.length && lanes[l] > leftEdge - 3) l++;
          lanes[l] = v._x + w / 2;
          v._lane = l;
        } else {
          v._lane = -1;
        }
      });
      // apply manual overrides (horizontal nudge dx, optional lane bump)
      list.forEach(function (v) {
        var o = ovFor(v);
        v._dx = (o && o.dx) || 0;
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
      var topAboveBox = W + maxLen + (maxLane >= 0 ? LABEL_H : 0) + PAD;
      sr.container.style.height = Math.max(0, topAboveBox - gap) + "px";
      list.forEach(function (v) { draw(sr, v); });
    });

    buildToggles();
    setupPosTip(table);
    window.__variantInfo = { total: DATA.length, missing: missing, counts: counters };
    if (missing) console.warn("variant renderer: " + missing + " unmapped variants");
  }

  function draw(sr, v) {
    var top = (sr.side === "top");
    var gap = sr._gap || 0;
    var W = whiteGap();              // white space before the AA box
    var len = TICK + (v._lane >= 0 ? v._lane * LANE_H : 0);  // FIXED line length
    var base = gap - W;              // line bottom sits W above the actual letter box

    var tick = document.createElement("div");
    tick.className = "vtick";
    tick.style.left = (v._x - 1) + "px";
    tick.style.height = len + "px";
    tick.style.background = v._color;
    tick.style.color = v._color;   // for the hover glow (currentColor)
    tick.style[top ? "bottom" : "top"] = (-base) + "px";
    setData(tick, v);
    sr.container.appendChild(tick);

    if (v._lane >= 0) {
      var lab = document.createElement("div");
      lab.className = "vlabel";
      lab.textContent = v.label;
      lab.style.left = (v._x + v._dx) + "px";
      lab.style.color = v._color;
      lab.style[top ? "bottom" : "top"] = (-base + len + 1) + "px";
      setData(lab, v);
      sr.container.appendChild(lab);
    }
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
      else if (!(e.target.closest && e.target.closest("#variant-popup"))) hide(pop);
    });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") hide(pop); });
    window.addEventListener("resize", function () { hide(pop); });
  }
  function hide(pop) { pop.style.display = "none"; }

  function openPopup(marker, pop) {
    var protein = marker.getAttribute("data-protein");
    var mut = marker.getAttribute("data-mutation");
    var cat = marker.getAttribute("data-category");
    var origin = marker.getAttribute("data-origin");
    var color = (CFG[cat] || CFG.Other).color;
    var d = DETAILS[protein + ":" + mut] || DETAILS[mut];

    var h = '<button class="vpop-close" aria-label="Close">&times;</button>';
    h += '<div class="vpop-head"><span class="vpop-chip" style="background:' + color +
         '"></span>' + esc(mut) + '<span class="vpop-prot">' + esc(protein) + '</span></div>';

    if (d) {
      var paper = (d.paper && PAPERS[d.paper]) || {};
      var url = paper.url || (paper.pmid ? "https://pubmed.ncbi.nlm.nih.gov/" + paper.pmid + "/" : null);
      if (paper.title) {
        h += '<a class="vpop-title" ' + (url ? 'href="' + esc(url) + '" target="_blank" rel="noopener"' : "") +
             ' title="' + esc(paper.title) + '">' + esc(paper.title) + "</a>";
      }
      if (paper.pmid) {
        h += '<div class="vpop-row vpop-dim">PMID ' +
             '<a href="https://pubmed.ncbi.nlm.nih.gov/' + esc(paper.pmid) + '/" target="_blank" rel="noopener">' +
             esc(paper.pmid) + "</a></div>";
      }
      if (d.gnomad) {
        var g = d.gnomad.present
          ? "Yes" + (d.gnomad.maf ? " (MAF " + esc(d.gnomad.maf) + ")" : "")
          : "No";
        h += '<div class="vpop-row"><b>gnomAD:</b> ' + g + "</div>";
      }
      var eff = d.effect || cat;
      h += '<div class="vpop-row"><b>' + esc(eff) + "</b>" +
           (d.method ? ' <span class="vpop-dim">(' + esc(d.method) + ")</span>" : "") + "</div>";
    } else {
      h += '<div class="vpop-row vpop-dim">Details not added yet' +
           (origin ? " &middot; source: " + esc(origin) : "") + "</div>";
    }
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
    });
    document.body.appendChild(box);
  }

  var t;
  function schedule() { clearTimeout(t); t = setTimeout(run, 120); }

  // Load the editable data files fresh each time so "edit a data_*.js + refresh"
  // always shows your changes (browsers otherwise cache <script src>). On a
  // double-clicked file:// page we skip the query (not needed, and keeps it simple).
  function boot() {
    var files = ["data_variants.js", "data_overrides.js", "data_papers.js", "data_details.js"];
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
