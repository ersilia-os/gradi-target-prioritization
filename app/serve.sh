#!/usr/bin/env bash
# Serve the target-prioritization browser on the LAN/Tailscale.
# Bound to 0.0.0.0 so it is reachable from other tailnet devices (e.g. your laptop):
#   http://miquel-macmini:8080   or   http://100.105.44.112:8080
# Locally on this Mac mini: http://localhost:8080
#
# Uses the gradi env python (stdlib http.server; no extra deps). Override PORT via env.
set -euo pipefail
cd "$(dirname "$0")"

PORT="${PORT:-8080}"
PY="${GRADI_PY:-$HOME/miniconda3/envs/gradi/bin/python}"
[ -x "$PY" ] || PY="python3"

echo "Serving $(pwd) on http://0.0.0.0:${PORT}"
echo "  local:     http://localhost:${PORT}"
echo "  tailscale: http://miquel-macmini:${PORT}  (or http://100.105.44.112:${PORT})"
exec "$PY" -m http.server "$PORT" --bind 0.0.0.0
