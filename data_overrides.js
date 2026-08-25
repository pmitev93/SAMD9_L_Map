/* Manual label nudges, keyed by label (or "PROTEIN:label").
   dx    = pixels right(+)/left(-)
   dlen  = extra line length in px (label moves out with the line; + = longer)
   lane  = stack level (0 = shortest line)
   color = override this one variant's tick/label color (any CSS color) */
window.VARIANT_OVERRIDES = {
  "Y72C": { "lane": 0, "dx": 14 },
  "R70C": { "dx": -14 },
  "A1195V": { "dx": -14 },
  "V1276I": { "dlen": 10 },
  "H464Y": { "color": "orange" },
  "H468Y": { "color": "orange" }
};
