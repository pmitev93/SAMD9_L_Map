/* Data-driven variant renderer for the SAMD9/SAMD9L map.
   Reads embedded variants.json, removes the old hand-placed markers, and draws
   colored ticks + lane-stacked labels above SAMD9L / below SAMD9. */
(function () {
  "use strict";

  var CFG = {
    GoF:    { color: "#e02424", label: true,  legend: "GoF (gain of function)" },
    LoF:    { color: "#1f56bc", label: true,  legend: "LoF (loss of function)" },
    gnomAD: { color: "#3aa83a", label: false, legend: "gnomAD (truncating)" },
    Other:  { color: "#333333", label: true,  legend: "Other (somatic, NoF, ...)" }
  };
  var TICK = 16, LANE_H = 13, LABEL_H = 13, PAD = 4;

  function estW(s) { return s.length * 6.1 + 6; }

  function run() {
    var dataEl = document.getElementById("variant-data");
    var table = document.querySelector("table");
    if (!dataEl || !table) return;
    var DATA = JSON.parse(dataEl.textContent);

    // 1. remove old hand-placed markers
    table.querySelectorAll(".vertical-line").forEach(function (n) { n.remove(); });
    // clear any previous render (idempotent on resize)
    table.querySelectorAll("tr.vmark").forEach(function (n) { n.remove(); });

    // 2. locate sequence rows and insert a reserved-height container next to each
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
      td.appendChild(div);
      tr.appendChild(td);
      if (side === "top") row.parentNode.insertBefore(tr, row);
      else row.parentNode.insertBefore(tr, row.nextSibling);
      seqRows.push({ protein: t, row: row, side: side, container: div });
    });

    // 3. residue -> cell map (cumulative per protein, in document order)
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
        }
      }
    });

    // 4. attach variants to their block group
    var groups = new Map();
    var missing = 0;
    DATA.forEach(function (v) {
      var cell = maps[v.protein] && maps[v.protein][v.residue];
      if (!cell) { missing++; return; }
      v._cell = cell;
      var sr = cellRow.get(cell);
      if (!groups.has(sr)) groups.set(sr, []);
      groups.get(sr).push(v);
    });

    // 5. lay out + draw each group
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
          if (l > maxLane) maxLane = l;
        } else {
          v._lane = -1;
        }
      });

      var stackH = (maxLane + 1) * LANE_H + (maxLane >= 0 ? LABEL_H : 0) + PAD;
      var H = TICK + stackH;
      sr.container.style.height = H + "px";

      list.forEach(function (v) { draw(sr, v); });
    });

    // 6. expose for toggles/tooltips (Stage 3)
    window.__variantInfo = { total: DATA.length, missing: missing,
                             counts: counters };
    if (missing) console.warn("variant renderer: " + missing + " unmapped variants");
  }

  function draw(sr, v) {
    var top = (sr.side === "top");
    var tickLen = TICK + (v._lane >= 0 ? v._lane * LANE_H : 0);

    var tick = document.createElement("div");
    tick.className = "vtick";
    tick.style.left = (v._x - 1) + "px";
    tick.style.height = tickLen + "px";
    tick.style.background = v._color;
    tick.style[top ? "bottom" : "top"] = "0px";
    setData(tick, v);
    sr.container.appendChild(tick);

    if (v._lane >= 0) {
      var lab = document.createElement("div");
      lab.className = "vlabel";
      lab.textContent = v.label;
      lab.style.left = v._x + "px";
      lab.style.color = v._color;
      lab.style[top ? "bottom" : "top"] = (tickLen + 1) + "px";
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
    el.title = v.label + "  —  " + v.protein + " · " +
               (v.effect || v.category) + (v.origin ? " · " + v.origin : "");
  }

  var t;
  function schedule() { clearTimeout(t); t = setTimeout(run, 120); }
  if (document.readyState === "loading")
    document.addEventListener("DOMContentLoaded", run);
  else run();
  window.addEventListener("resize", schedule);
})();
