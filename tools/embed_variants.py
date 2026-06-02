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
CSS  = ROOT / "tools" / "variant_styles.css"
JS   = ROOT / "tools" / "variant_renderer.js"

START = "<!-- VARIANT-LAYER-START -->"
END   = "<!-- VARIANT-LAYER-END -->"

def main():
    html = HTML.read_text(encoding="utf-8")
    data = DATA.read_text(encoding="utf-8")
    css  = CSS.read_text(encoding="utf-8")
    js   = JS.read_text(encoding="utf-8")

    block = (
        START + "\n"
        '<style>\n' + css + '\n</style>\n'
        '<script id="variant-data" type="application/json">\n' + data + '\n</script>\n'
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
