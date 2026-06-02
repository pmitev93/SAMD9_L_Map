/* Data-driven variant renderer for the SAMD9/SAMD9L map.
   Reads embedded variant-data (+ optional variant-overrides), removes the old
   hand-placed markers, and draws colored ticks + lane-stacked labels above
   SAMD9L / below SAMD9. Also builds the show/hide category toggles. */
(function () {
  "use strict";

  var CFG = {
    GoF:    { color: "#FF0000", label: true,  on: true,  legend: "GoF (gain of function)" },
    LoF:    { color: "#1f56bc", label: true,  on: true,  legend: "LoF (loss of function)" },
    gnomAD: { color: "#73d73c", label: false, on: false, legend: "gnomAD (truncating)" },
    Other:  { color: "#000000", label: true,  on: true,  legend: "Other (somatic / NoF)" }
  };
  // Geometry. LANE_H/LABEL_H track the label size so bigger labels still fit.
  var TICK = 16, LANE_H = 16, LABEL_H = 15, PAD = 5;

  function estW(s) { return s.length * 7.0 + 8; }   // ~ label pixel width

  function hasDigit(el) { return el && /\d/.test(el.textContent || ""); }

  function run() {
    var dataEl = document.getElementById("variant-data");
    var table = document.querySelector("table");
    if (!dataEl || !table) return;
    var DATA = JSON.parse(dataEl.textContent);
    var OV = {};
    var ovEl = document.getElementById("variant-overrides");
    if (ovEl) { try { OV = JSON.parse(ovEl.textContent || "{}"); } catch (e) {} }
    function ovFor(v) { return OV[v.label] || OV[v.protein + ":" + v.label] || null; }

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
        var prev = row.previousElementSibling;
        var anchor = hasDigit(prev) ? prev : row;       // the ruler, if present
        anchor.parentNode.insertBefore(tr, anchor);
      } else {
        var next = row.nextElementSibling;
        var a2 = hasDigit(next) ? next : row;
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

      var stackH = (maxLane + 1) * LANE_H + (maxLane >= 0 ? LABEL_H : 0) + PAD;
      sr.container.style.height = (TICK + stackH) + "px";
      // distance from the container edge to the sequence letters (the number
      // ruler sits in this gap); ticks extend across it so they touch the box.
      var cr = sr.container.getBoundingClientRect(), rr = sr.row.getBoundingClientRect();
      sr._gap = Math.max(0, Math.round(sr.side === "top" ? rr.top - cr.bottom
                                                         : cr.top - rr.bottom));
      list.forEach(function (v) { draw(sr, v); });
    });

    buildToggles();
    window.__variantInfo = { total: DATA.length, missing: missing, counts: counters };
    if (missing) console.warn("variant renderer: " + missing + " unmapped variants");
  }

  function draw(sr, v) {
    var top = (sr.side === "top");
    var gap = sr._gap || 0;
    var tickLen = TICK + (v._lane >= 0 ? v._lane * LANE_H : 0);

    var tick = document.createElement("div");
    tick.className = "vtick";
    tick.style.left = (v._x - 1) + "px";
    tick.style.height = (tickLen + gap) + "px";   // +gap reaches the letter box
    tick.style.background = v._color;
    tick.style[top ? "bottom" : "top"] = (-gap) + "px";
    setData(tick, v);
    sr.container.appendChild(tick);

    if (v._lane >= 0) {
      var lab = document.createElement("div");
      lab.className = "vlabel";
      lab.textContent = v.label;
      lab.style.left = (v._x + v._dx) + "px";
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
              '<span class="vt-sw" style="background:' + CFG[cat].color + '"></span>' +
              CFG[cat].legend + '</label>';
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
  if (document.readyState === "loading")
    document.addEventListener("DOMContentLoaded", run);
  else run();
  window.addEventListener("resize", schedule);
})();
