#!/bin/bash
# Multi-mutant Cartesian ddG: 8개 인간화 mutation 동시 적용
#
# WT 구조에 8개 FR mutation을 한번에 적용하여 전체 인간화의 ddG를 계산.
# 개별 mutation ddG의 합과 비교하여 epistatic effect 확인 가능.
#
# Input:
#   - ../01_structure_prediction/fold_anti_fap_nb_model_0.cif : WT AF3 structure
#   - 01_multi_mutant_ddg.txt : 8개 mutation (Rosetta format)
# Output:
#   - 01_multi_mutant.ddg : raw energy output
#
# Usage:
#   bash 01_run_multi_mutant_ddg.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STRUCTURE="${SCRIPT_DIR}/../01_structure_prediction/fold_anti_fap_nb_model_0.cif"
MUT_FILE="${SCRIPT_DIR}/01_multi_mutant_ddg.txt"

ROSETTA_BIN="/home/laugh/rosetta/source/bin/cartesian_ddg.linuxgccrelease"

if [ ! -f "$STRUCTURE" ]; then
    echo "ERROR: Structure not found: $STRUCTURE"
    exit 1
fi

echo "========================================"
echo "Multi-Mutant Cartesian ddG"
echo "========================================"
echo "Structure: $STRUCTURE"
echo "Mutations: 8 simultaneous (humanization)"
echo ""

cd "$SCRIPT_DIR"

"$ROSETTA_BIN" \
    -in:file:s "$STRUCTURE" \
    -ddg:mut_file "$MUT_FILE" \
    -ddg:iterations 5 \
    -ddg:cartesian \
    -ddg:bbnbrs 1 \
    -ddg:dump_pdbs false \
    -score:weights ref2015_cart \
    -fa_max_dis 9.0 \
    -ddg:legacy true \
    -out:file:scorefile 01_multi_mutant.sc \
    > 01_multi_mutant_run.log 2>&1

if [ -f "mutations.ddg" ]; then
    mv mutations.ddg 01_multi_mutant.ddg
    echo "Done. Output: 01_multi_mutant.ddg"
else
    echo "FAIL - check 01_multi_mutant_run.log"
fi
