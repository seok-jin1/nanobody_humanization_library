#!/bin/bash
# WT vs Humanized 구조 FastRelax 후 total score 비교
#
# 양쪽 AF3 구조를 Rosetta FastRelax로 에너지 최소화 한 뒤
# total_score를 비교하여 전반적 안정성을 평가한다.
#
# Input:
#   - ../01_structure_prediction/fold_anti_fap_nb_model_0.cif      : WT structure
#   - ../01_structure_prediction/fold_humanized_fap_nb_model_0.cif : Humanized structure
# Output:
#   - 02_relax_wt.sc          : WT relaxed score
#   - 02_relax_humanized.sc   : Humanized relaxed score
#
# Usage:
#   bash 02_run_relax_compare.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WT_STRUCTURE="${SCRIPT_DIR}/../01_structure_prediction/fold_anti_fap_nb_model_0.cif"
HUM_STRUCTURE="${SCRIPT_DIR}/../01_structure_prediction/fold_humanized_fap_nb_model_0.cif"

RELAX_BIN="/home/laugh/rosetta/source/bin/relax.linuxgccrelease"

if [ ! -f "$WT_STRUCTURE" ]; then
    echo "ERROR: WT structure not found: $WT_STRUCTURE"
    exit 1
fi
if [ ! -f "$HUM_STRUCTURE" ]; then
    echo "ERROR: Humanized structure not found: $HUM_STRUCTURE"
    exit 1
fi

echo "========================================"
echo "FastRelax: WT vs Humanized"
echo "========================================"

cd "$SCRIPT_DIR"

# WT relax
echo "[$(date +%H:%M:%S)] Relaxing WT structure..."
"$RELAX_BIN" \
    -in:file:s "$WT_STRUCTURE" \
    -relax:cartesian \
    -score:weights ref2015_cart \
    -nstruct 5 \
    -out:file:scorefile 02_relax_wt.sc \
    -out:prefix wt_ \
    > 02_relax_wt.log 2>&1
echo "[$(date +%H:%M:%S)] WT relax done"

# Humanized relax
echo "[$(date +%H:%M:%S)] Relaxing Humanized structure..."
"$RELAX_BIN" \
    -in:file:s "$HUM_STRUCTURE" \
    -relax:cartesian \
    -score:weights ref2015_cart \
    -nstruct 5 \
    -out:file:scorefile 02_relax_humanized.sc \
    -out:prefix hum_ \
    > 02_relax_humanized.log 2>&1
echo "[$(date +%H:%M:%S)] Humanized relax done"

echo ""
echo "========================================"
echo "Results:"
echo "  WT scores:        02_relax_wt.sc"
echo "  Humanized scores: 02_relax_humanized.sc"
echo "========================================"
echo ""
echo "Next: python 03_analyze_relax.py"
