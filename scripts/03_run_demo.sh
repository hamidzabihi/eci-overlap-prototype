#!/usr/bin/env bash
set -euo pipefail
: "${ECI_PROTO:?Set ECI_PROTO}"
exec "$ECI_PROTO/examples/BA_BD_overlap/run.sh" "$@"
