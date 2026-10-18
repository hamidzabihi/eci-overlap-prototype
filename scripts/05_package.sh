#!/usr/bin/env bash
set -euo pipefail
: "${ECI_PROTO:?Set ECI_PROTO}"

OUT="$ECI_PROTO/../eci2game_phase1_$(date +%Y%m%d).tar.gz"
cd "$ECI_PROTO"

for f in bin/SHARC_ECI_GAME.py \
         examples/BA_BD_overlap/run.sh \
         examples/BA_BD_overlap/QM.in.step0 \
         examples/BA_BD_overlap/QM.in.step1 \
         docs/ECI2GAME_status.md \
         docs/upstream_pr.md docs/NOTES_qmout_validation.md scripts/read_overlap.py ; do
    [[ -e "$f" ]] || { echo "missing $f - run earlier scripts first" >&2; exit 1; }
done

if [[ ! -f examples/BA_BD_overlap/QM.out.step1 ]]; then
    echo "  warn: QM.out.step1 missing (no QC backend on this box)"
    echo "        tarball will not contain reference outputs"
    echo "        Tomislav will regenerate them by running run.sh"
fi

tar czf "$OUT" \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.pytest_cache' \
    --exclude='examples/BA_BD_overlap/SAVE' \
    --exclude='examples/BA_BD_overlap/*.out.*' \
    src/ scripts/ bin/ examples/BA_BD_overlap/ tests/ docs/ \
    PHASE1.md README.md NOTES.md NEXT_MILESTONE.md CITATION.cff pytest.ini \
    2>/dev/null || \
tar czf "$OUT" \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.pytest_cache' \
    src/ scripts/ bin/ examples/BA_BD_overlap/ tests/ docs/ \
    README.md CITATION.cff pytest.ini

echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
echo "== top of tarball =="
tar tzf "$OUT" | head -20
