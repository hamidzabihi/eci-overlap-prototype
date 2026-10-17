#!/usr/bin/env bash
set -euo pipefail
: "${SHARC:?Set SHARC to the sharc4 install directory}"
: "${ECI_PROTO:?Set ECI_PROTO to the eci-overlap-prototype directory}"

die() { echo "FAIL: $*" >&2; exit 1; }
ok()  { echo "  ok: $*"; }

echo "== environment check =="

python -c 'import sys; assert sys.version_info >= (3,10)' \
    || die "Python 3.10+ required"
ok "Python $(python -c 'import sys;print(sys.version.split()[0])')"

for pkg in numpy scipy yaml pyscf; do
    python -c "import $pkg" 2>/dev/null || die "missing python package: $pkg"
    ok "python package: $pkg"
done

[[ -f "$SHARC/bin/SHARC_ECI.py" ]] || die "SHARC_ECI.py not found in $SHARC/bin"
ok "SHARC_ECI.py at $SHARC/bin"

[[ -f "$ECI_PROTO/src/esd_overlap.py" ]] || die "missing esd_overlap.py"
[[ -f "$ECI_PROTO/src/state_overlap.py" ]] || die "missing state_overlap.py"
[[ -f "$ECI_PROTO/src/phase_tracking.py" ]] || die "missing phase_tracking.py"
ok "prototype modules present"

export PYTHONPATH="$SHARC/bin:$SHARC/lib:${PYTHONPATH:-}"
python -c 'import SHARC_INTERFACE' || die "cannot import SHARC_INTERFACE"
python -c 'import SHARC_ECI'       || die "cannot import SHARC_ECI"
ok "sharc4 libraries importable"

if command -v g16 >/dev/null 2>&1 || command -v gaussian >/dev/null 2>&1; then
    ok "Gaussian found"
elif command -v orca >/dev/null 2>&1; then
    ok "ORCA found"
else
    echo "  warn: no Gaussian or ORCA on PATH; demo run will fail"
fi

echo "== environment ok =="
