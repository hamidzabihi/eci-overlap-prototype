#!/usr/bin/env bash
set -euo pipefail
: "${ECI_PROTO:?Set ECI_PROTO}"

DST="$ECI_PROTO/examples/BA_BD_overlap"
SKIP=0
[[ -f "$DST/QM.out.step0" ]] || SKIP=1
[[ -f "$DST/QM.out.step1" ]] || SKIP=1

if [[ $SKIP -eq 1 ]]; then
    echo "== verify: skipped =="
    echo "  QM.out.step0 and/or QM.out.step1 not yet produced."
    echo "  Run 03_run_demo.sh on a machine with a working QC backend"
    echo "  (Gaussian or ORCA, per ECI.template). Format and scaffolding"
    echo "  are ready; only the numerical verification is deferred."
    exit 0
fi

fail=0
note()  { echo "  $*"; }
check() { if eval "$1"; then note "ok: $2"; else note "FAIL: $2"; fail=1; fi; }

echo "== step 0 =="
check "! grep -q '^! 6 ' '$DST/QM.out.step0'" "no ! 6 section (init step)"
check "! grep -q '^! 7 ' '$DST/QM.out.step0'" "no ! 7 section (init step)"
check "grep -q '^! 1 ' '$DST/QM.out.step0'"  "has ! 1 Hamiltonian"
check "grep -q '^! 8 ' '$DST/QM.out.step0'"  "has ! 8 runtime"

echo "== step 1 =="
check "grep -q '^! 6 Overlap matrix' '$DST/QM.out.step1'"       "has ! 6 Overlap"
check "grep -q '^! 7 Wave function phases' '$DST/QM.out.step1'" "has ! 7 phases"

if [[ $fail -eq 0 ]]; then
    python - "$DST/QM.out.step1" <<'PY'
import sys, re, numpy as np
path = sys.argv[1]
def parse(text, flag):
    m = re.search(rf"! {flag}[^\n]*\n(\d+) (\d+)\n(.*?)\n\n", text, re.S)
    if not m:
        raise SystemExit(f"! {flag} not parseable")
    n = int(m.group(1)); vals = m.group(3).split()
    arr = np.array([complex(float(vals[2*i]), float(vals[2*i+1]))
                    for i in range(n*n)])
    return arr.reshape(n, n)
txt = open(path).read()
S = parse(txt, 6)
diag = np.diagonal(S)
print(f"  overlap shape:     {S.shape}")
print(f"  |diag| min/max:    {np.abs(diag).min():.6f} / {np.abs(diag).max():.6f}")
print(f"  max |imag(diag)|:  {np.abs(diag.imag).max():.2e}")
assert np.allclose(np.abs(diag), 1.0, atol=1e-6), "diag not unit modulus"
assert np.abs(diag.imag).max() < 1e-6, "diag has non-negligible Im"
if "! 7 " in txt:
    P = parse(txt, 7)
    assert np.allclose(np.abs(P.diagonal()), 1.0, atol=1e-12)
    print(f"  phases shape:      {P.shape}, all |phase|=1")
print("  numerical checks:  PASS")
PY
    rc=$?
    [[ $rc -ne 0 ]] && fail=1
fi

exit $fail
