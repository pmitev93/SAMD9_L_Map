#!/bin/bash
# Double-click to preview locally with live JSON (edit a .json + refresh = updates).
cd "$(dirname "$0")"
PORT=8765
URL="http://localhost:$PORT/Comprehensive_map_ver2.html"
echo "Serving at $URL   (Ctrl-C to stop)"
( sleep 1; open "$URL" ) &
python3 -m http.server $PORT
