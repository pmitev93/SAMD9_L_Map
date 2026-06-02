#!/bin/bash
# Double-click this after editing the Excel master to rebuild the map.
# (Edits the variant list from Conservation_Mutational_Landscape_Both.xlsx,
#  then embeds it into Comprehensive_map_ver2.html.)
cd "$(dirname "$0")"
echo "Rebuilding variants from the Excel master..."
python3 tools/build_variants.py && python3 tools/embed_variants.py
echo ""
echo "Done. Reload Comprehensive_map_ver2.html in your browser."
read -p "Press Enter to close."
