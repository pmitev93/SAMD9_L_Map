#!/usr/bin/env python3
"""Embed variants.json + the variant CSS/JS into Comprehensive_map_ver2.html.

Idempotent: replaces the block between the VARIANT-LAYER markers if present,
otherwise injects it just before the page's trailing <script>.

Run:  python3 tools/embed_variants.py   (run build_variants.py first)
"""
import pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
HTML = ROOT / "Comprehensive_map_ver2.html"
DATA = ROOT / "variants.json"
OVER = ROOT / "variant_overrides.json"
CSS  = ROOT / "tools" / "variant_styles.css"
JS   = ROOT / "tools" / "variant_renderer.js"

# Manual label nudges survive re-running build_variants.py because they live here,
# not in variants.json. Key by label (or "PROTEIN:label"); dx = pixels right(+)/
# left(-); lane optionally bumps a label up/down a level. Example:
#   { "S1550P": { "dx": -14 }, "SAMD9:R902W": { "dx": 10, "lane": 1 } }
DEFAULT_OVERRIDES = "{\n}\n"

START = "<!-- VARIANT-LAYER-START -->"
END   = "<!-- VARIANT-LAYER-END -->"

def main():
    html = HTML.read_text(encoding="utf-8")
    data = DATA.read_text(encoding="utf-8")
    css  = CSS.read_text(encoding="utf-8")
    js   = JS.read_text(encoding="utf-8")
    if not OVER.exists():
        OVER.write_text(DEFAULT_OVERRIDES, encoding="utf-8")
    over = OVER.read_text(encoding="utf-8")

    block = (
        START + "\n"
        '<style>\n' + css + '\n</style>\n'
        '<script id="variant-data" type="application/json">\n' + data + '\n</script>\n'
        '<script id="variant-overrides" type="application/json">\n' + over + '\n</script>\n'
        '<script>\n' + js + '\n</script>\n'
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
    print(f"Variant layer {action} ({len(data)} bytes data).")

if __name__ == "__main__":
    main()
