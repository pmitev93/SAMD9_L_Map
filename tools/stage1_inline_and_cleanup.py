#!/usr/bin/env python3
"""Stage 1: make Comprehensive_map_ver2.html self-contained + cleaner.

Reproducible transform (reads asset files, edits the HTML in place):
  1. Inline Mitev_MSA_Colored.NEW.EM.css into a <style> block (drop the bogus
     trailing <style>/comment junk in that file). Remove the <link>.
  2. Append "enhancement" CSS: modern system font stack, a real spacer row to
     replace the white-"B" hack, and styling for the rebuilt legends.
  3. Replace the white-"B" spacer cells (font-size:5px/2px) with empty .seq-gap
     cells.
  4. Replace the bottom legend: rebuild the ConSurf scale in HTML/CSS, inline the
     Mutations legend SVG, keep the Domains scheme as an external (swappable) file.

Run from the project folder:  python3 tools/stage1_inline_and_cleanup.py
"""
import sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
HTML = ROOT / "Comprehensive_map_ver2.html"
CSS  = ROOT / "Mitev_MSA_Colored.NEW.EM.css"
MUT  = ROOT / "Mutations_Legend.svg"

def expect(cond, msg):
    if not cond:
        sys.exit(f"ABORT: {msg}")

html = HTML.read_text(encoding="utf-8")

# --- 1. clean the external CSS (keep everything before the stray <style>) ------
css = CSS.read_text(encoding="utf-8")
expect("\n<style>" in css, "expected stray <style> marker in CSS file")
base_css = css.split("\n<style>", 1)[0].rstrip()

# --- 2. enhancement CSS --------------------------------------------------------
ENHANCE = r"""
/* ===================== Stage 1 enhancements ===================== */
:root{
  --ui-font: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue",
             Arial, sans-serif;
}
body{ font-family: var(--ui-font); }
h1{ font-family: var(--ui-font); font-weight: 600; letter-spacing: .01em; }
/* modern font everywhere, overriding legacy <font face="..."> attributes.
   Column alignment is preserved because the table is table-layout:fixed with
   fixed 1em-wide, centered cells (alignment comes from cell width, not glyph). */
td, td font, font{ font-family: var(--ui-font); }

/* Proper vertical gap between the SAMD9L and SAMD9 rows.
   Replaces the old white-on-white "B" character hack. */
tr.seq-gap td, td.seq-gap{ height: 6px; line-height: 6px; padding: 0; font-size: 0; }

/* ---- Legends (bottom bar) ---- */
.legend-block{ display:flex; flex-direction:column; align-items:center;
               margin:0 28px; font-family:var(--ui-font); }
.consurf-row{ display:flex; }
.cs-box{ width:22px; height:22px; box-sizing:border-box; display:inline-flex;
         align-items:center; justify-content:center; font-size:13px;
         font-weight:600; border:1px solid #cccccc; }
.consurf-labels{ display:flex; justify-content:space-between; width:198px;
                 font-size:11px; color:#333; margin-top:3px; }
.consurf-isd{ display:flex; align-items:center; gap:6px; font-size:11px;
              color:#333; margin-top:5px; }
.cs-isd{ width:18px; height:14px; background:#ffff99; border:1px solid #cccccc;
         display:inline-block; }
"""

# --- inline CSS at the <link> location, append enhancements --------------------
LINK = '    <link rel="stylesheet" type="text/css" href="Mitev_MSA_Colored.NEW.EM.css">'
expect(html.count(LINK) == 1, "CSS <link> line not found exactly once")
html = html.replace(
    LINK,
    "    <style>\n" + base_css + "\n" + ENHANCE + "\n    </style>",
)

# --- 3. replace the white-"B" spacer cells ------------------------------------
b5 = '<td class="Score11" style="font-size: 5px;">B</td>'
b2 = '<td class="Score11" style="font-size: 2px;">B</td>'
n5, n2 = html.count(b5), html.count(b2)
expect(n5 == 31 and n2 == 1, f"unexpected B-spacer counts: 5px={n5}, 2px={n2}")
html = html.replace(b5, '<td class="seq-gap"></td>').replace(
    b2, '<td class="seq-gap"></td>')

# --- 4. rebuild the bottom legends --------------------------------------------
# ConSurf scale: reuse the existing Score1..Score9 background colors.
cs_cells = "".join(
    f'<span class="cs-box Score{i}">{i}</span>' for i in range(1, 10))
CONSURF_HTML = (
    '<div class="legend-block">'
    f'<div class="consurf-row">{cs_cells}</div>'
    '<div class="consurf-labels"><span>Variable</span><span>Average</span>'
    '<span>Conserved</span></div>'
    '<div class="consurf-isd"><span class="cs-isd"></span>Insufficient Data</div>'
    '</div>'
)

img_consurf = ('<img src="ConSurf_Color_Scale.png" alt="Italian Trulli" '
               'style="width:200px;height:50px; margin-left:10px;">')
expect(html.count(img_consurf) == 1, "ConSurf <img> not found exactly once")
html = html.replace(img_consurf, CONSURF_HTML)

# Mutations legend: inline the SVG (add viewBox + display size).
mut = MUT.read_text(encoding="utf-8").strip()
old_svg_open = '<svg width="1888" height="578" '
expect(old_svg_open in mut, "mutations SVG header not as expected")
mut = mut.replace(
    old_svg_open,
    '<svg viewBox="0 0 1888 578" style="height:80px;width:auto;display:block;" ',
    1,
)
img_mut = ('<img src="Mutations_Legend.svg" alt="Italian Trulli" '
           'style="height:80px; margin-left:40px;"> ')
expect(html.count(img_mut) == 1, "Mutations <img> not found exactly once")
html = html.replace(img_mut, '<div class="legend-block">' + mut + '</div>')

# Domains scheme: keep external (user will swap the file), just wrap for spacing.
img_dom = ('<img src="Domains_Scheme.svg" alt="Italian Trulli" '
           'style="height:80px; margin-left:40px;">')
expect(html.count(img_dom) == 1, "Domains <img> not found exactly once")
html = html.replace(
    img_dom,
    '<div class="legend-block"><img src="Domains_Scheme.svg" '
    'alt="Domain scheme" style="height:80px;"></div>',
)

HTML.write_text(html, encoding="utf-8")
print("Stage 1 transforms applied OK.")
print(f"  inlined CSS ({len(base_css)} chars) + enhancements")
print(f"  replaced {n5+n2} white-B spacer cells")
print("  rebuilt ConSurf scale (HTML), inlined Mutations SVG, kept Domains as file")
