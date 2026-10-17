#!/usr/bin/env bash
set -euo pipefail
: "${SHARC:?Set SHARC}" "${ECI_PROTO:?Set ECI_PROTO}"

SRC="$SHARC/examples/SHARC_ECI"
DST="$ECI_PROTO/examples/BA_BD_overlap"
[[ -d "$SRC" ]] || { echo "missing $SRC" >&2; exit 1; }
[[ -f "$SRC/QM.in" ]] || { echo "missing $SRC/QM.in" >&2; exit 1; }

rm -rf "$DST"
mkdir -p "$DST"

for f in ECI.template ECI.resources QM.in; do
    [[ -f "$SRC/$f" ]] && cp "$SRC/$f" "$DST/"
done

NATOM=$(head -1 "$SRC/QM.in")
GEOM_LINES=$((NATOM + 2))
head -n "$GEOM_LINES" "$SRC/QM.in" > "$DST/_geom.xyz"

# Step 0: initial, no overlap (SHARC forbids overlap at step 0)
cp "$DST/_geom.xyz" "$DST/QM.in.step0"
cat >> "$DST/QM.in.step0" <<'EOF'
step 0
unit angstrom
states 6 0 6
charge 0 1 0
h
dm
EOF

# Step 1: same atoms, tiny displacement on atom 1, with overlap + phases
cp "$DST/_geom.xyz" "$DST/QM.in.step1"
python - "$DST/QM.in.step1" <<'PY'
import sys
lines = open(sys.argv[1]).read().splitlines()
parts = lines[2].split()
parts[1] = f"{float(parts[1]) + 0.01:.6f}"
lines[2] = "  ".join(parts)
open(sys.argv[1], "w").write("\n".join(lines) + "\n")
PY
cat >> "$DST/QM.in.step1" <<'EOF'
step 1
unit angstrom
states 6 0 6
charge 0 1 0
h
dm
overlap
phases
EOF

cat > "$DST/run.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
: "${SHARC:?Set SHARC}" "${ECI_PROTO:?Set ECI_PROTO}"
export PYTHONPATH="$SHARC/bin:$SHARC/lib:$ECI_PROTO/bin:${PYTHONPATH:-}"
rm -rf SAVE

echo "== step 0 (h + dm) =="
cp QM.in.step0 QM.in
python -m SHARC_ECI_GAME QM.in 2>&1 | tail -n 30
mv QM.out QM.out.step0

echo "== step 1 (h + dm + overlap + phases) =="
cp QM.in.step1 QM.in
python -m SHARC_ECI_GAME QM.in 2>&1 | tail -n 30
mv QM.out QM.out.step1

echo "== done: QM.out.step0, QM.out.step1 =="
EOF
chmod +x "$DST/run.sh"

cat > "$DST/README.md" <<'MDEOF'
# BA/BD overlap demo

Two-step run of `SHARC_ECI_GAME` on the BA/BD two-fragment example
distributed with `sharc4/examples/SHARC_ECI`.

- `QM.in.step0` - h + dm, no overlap (SHARC forbids overlap at step 0)
- `QM.in.step1` - h + dm + overlap + phases; geometry shifted by 0.01 A
  on the first atom so the overlap matrix is non-trivial

## Requirements

- A working sharc4 install: `export SHARC=/path/to/sharc4`
- This repo's `bin/` on PYTHONPATH: `export ECI_PROTO=$(pwd)/..`
- A working Gaussian installation (the example's `ECI.template` points
  at the GAUSSIAN interface). For an ORCA-based demo, edit
  `ECI.template` and change every `interface: GAUSSIAN` to
  `interface: ORCA`; also supply an `ORCA.resources` and rewrite
  `ECI.resources` accordingly.

## Run

    ./run.sh

Produces `QM.out.step0` (sections `! 1`, `! 2`, `! 8`) and
`QM.out.step1` (also `! 6 Overlap matrix` and `! 7 Wavefunction phases`).

## What to expect

- `! 6` is a 24x24 complex matrix, diagonal ≈ 1 (unit modulus, real positive).
- `! 7` is a 24-element vector of complex phase factors with |phase| = 1.

If you see `NotImplementedError: ECI basis changed for mult ...`, the
ECI basis expanded or contracted between steps 0 and 1. This is
Phase-2 scope, not a bug: the prototype only handles basis-stable
overlaps. Reduce the step-1 displacement or relax `tO` in
`ECI.template` until the basis stabilises.
MDEOF

rm -f "$DST/_geom.xyz"
echo "wrote $DST"
ls "$DST"
