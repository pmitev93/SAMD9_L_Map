#!/usr/bin/env python3
"""Embed the variant CSS + renderer (code) into Comprehensive_map_ver2.html, and
wire up the data files via <script src> tags.

The DATA itself lives in editable data_*.js files (window.VARIANTS / PAPERS /
VARIANT_DETAILS / VARIANT_OVERRIDES). Those load via <script src>, which works
even on a double-clicked file:// page — so you edit a data file, refresh, done
(no rebuild). You only re-run THIS script if the renderer code/CSS changes.

Idempotent: replaces the block between the VARIANT-LAYER markers if present,
otherwise injects it just before the page's trailing <script>.

Run:  python3 tools/embed_variants.py
"""
import pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
HTML = ROOT / "Comprehensive_map_ver2.html"
CSS  = ROOT / "tools" / "variant_styles.css"
JS   = ROOT / "tools" / "variant_renderer.js"

START = "<!-- VARIANT-LAYER-START -->"
END   = "<!-- VARIANT-LAYER-END -->"

def main():
    html = HTML.read_text(encoding="utf-8")
    css  = CSS.read_text(encoding="utf-8")
    js   = JS.read_text(encoding="utf-8")

    block = (
        START + "\n"
        "<style>\n" + css + "\n</style>\n"
        + "<script>\n" + js + "\n</script>\n"
        + END
    )

    if START in html and END in html:
        html = re.sub(re.escape(START) + r".*?" + re.escape(END), lambda m: block,
                      html, flags=re.S)
        action = "replaced"
    else:
        anchor = "<script>\n    window.onload = function() {"
        if anchor not in html:
            sys.exit("ABORT: injection anchor not found")
        html = html.replace(anchor, block + "\n" + anchor, 1)
        action = "injected"

    HTML.write_text(html, encoding="utf-8")
    print("Variant layer %s (renderer+CSS embedded; data via <script src>)." % action)

if __name__ == "__main__":
    main()
