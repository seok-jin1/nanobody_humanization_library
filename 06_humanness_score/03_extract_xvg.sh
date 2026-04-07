#!/bin/bash
# WT / Humanized MD trajectory에서 RMSD, RMSF, Rg xvg 파일 추출
#
# 03_compare_md.py 실행 전에 이 스크립트로 xvg 파일을 생성해야 함.
# WT trajectory 경로는 실제 위치에 맞게 수정 필요.
#
# Usage:
#   bash 03_extract_xvg.sh

set -euo pipefail

GMX=/usr/local/gromacs/bin/gmx
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# WT trajectory — 실제 경로로 수정
WT_TPR="/home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa/md.tpr"
WT_XTC="/home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa/md.xtc"
WT_OUT="${SCRIPT_DIR}/md_wt_xvg"

# Humanized trajectory
HUM_TPR="${SCRIPT_DIR}/md_humanized/md.tpr"
HUM_XTC="${SCRIPT_DIR}/md_humanized/md.xtc"
HUM_OUT="${SCRIPT_DIR}/md_humanized"

mkdir -p "$WT_OUT"

extract_analysis() {
    local tpr="$1"
    local xtc="$2"
    local outdir="$3"
    local label="$4"

    if [ ! -f "$tpr" ] || [ ! -f "$xtc" ]; then
        echo "[$label] TPR or XTC not found, skipping"
        return
    fi

    echo "[$label] Extracting RMSD..."
    echo -e "4\n4" | $GMX rms -s "$tpr" -f "$xtc" -o "$outdir/rmsd.xvg" -tu ns 2>/dev/null

    echo "[$label] Extracting RMSF..."
    echo "3" | $GMX rmsf -s "$tpr" -f "$xtc" -o "$outdir/rmsf.xvg" -res 2>/dev/null

    echo "[$label] Extracting Rg..."
    echo "1" | $GMX gyrate -s "$tpr" -f "$xtc" -o "$outdir/gyrate.xvg" 2>/dev/null

    echo "[$label] Done"
}

echo "========================================"
echo "Extracting MD analysis (xvg files)"
echo "========================================"

extract_analysis "$WT_TPR" "$WT_XTC" "$WT_OUT" "WT"
extract_analysis "$HUM_TPR" "$HUM_XTC" "$HUM_OUT" "Humanized"

echo ""
echo "Next: python 03_compare_md.py"
