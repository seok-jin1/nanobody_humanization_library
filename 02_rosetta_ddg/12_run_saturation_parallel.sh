#!/bin/bash
# FR saturation mutagenesis 병렬 실행 스크립트
#
# 각 FR position별 mutation file을 Rosetta cartesian_ddg로 병렬 실행.
# GNU parallel 또는 xargs를 사용하여 N_JOBS 만큼 동시 실행.
#
# Input:
#   - saturation/mutations_fr_saturation_split/*.txt : position별 mutation files
#   - ../01_structure_prediction/fold_anti_fap_nb_model_0.cif : AF3 structure
# Output:
#   - saturation/ddg_results/*.ddg : position별 ddG raw output
#
# Usage:
#   bash 12_run_saturation_parallel.sh [N_JOBS]
#   bash 12_run_saturation_parallel.sh 6     # 6코어 병렬

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STRUCTURE="${SCRIPT_DIR}/../01_structure_prediction/fold_anti_fap_nb_model_0.cif"
MUT_DIR="${SCRIPT_DIR}/saturation/mutations_fr_saturation_split"
OUT_DIR="${SCRIPT_DIR}/saturation/ddg_results"
LOG_DIR="${SCRIPT_DIR}/saturation/logs"

N_JOBS=${1:-4}  # default 4 parallel jobs

# Rosetta binary - adjust path as needed
ROSETTA_BIN="cartesian_ddg.linuxgccrelease"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Check prerequisites
if ! command -v "$ROSETTA_BIN" &>/dev/null; then
    echo "ERROR: $ROSETTA_BIN not found in PATH"
    echo "Set ROSETTA_BIN variable or add to PATH"
    exit 1
fi

if [ ! -f "$STRUCTURE" ]; then
    echo "ERROR: Structure file not found: $STRUCTURE"
    exit 1
fi

# Count mutation files
MUT_FILES=("$MUT_DIR"/*.txt)
N_FILES=${#MUT_FILES[@]}
echo "========================================"
echo "FR Saturation Mutagenesis - Parallel Run"
echo "========================================"
echo "Structure: $STRUCTURE"
echo "Mutation files: $N_FILES"
echo "Parallel jobs: $N_JOBS"
echo "Output dir: $OUT_DIR"
echo ""

# Function to run one mutation file
run_one() {
    local mut_file="$1"
    local basename=$(basename "$mut_file" .txt)
    local out_ddg="${OUT_DIR}/${basename}.ddg"
    local log_file="${LOG_DIR}/${basename}.log"

    # Skip if already completed
    if [ -f "$out_ddg" ]; then
        echo "[SKIP] $basename (already exists)"
        return 0
    fi

    echo "[START] $basename ($(date +%H:%M:%S))"

    # Run in a temp directory to avoid file conflicts
    local tmpdir=$(mktemp -d "/tmp/ddg_${basename}_XXXXX")
    cp "$mut_file" "$tmpdir/mutations.txt"

    cd "$tmpdir"
    "$ROSETTA_BIN" \
        -in:file:s "$STRUCTURE" \
        -ddg:mut_file mutations.txt \
        -ddg:iterations 5 \
        -ddg:cartesian \
        -ddg:bbnbrs 1 \
        -ddg:dump_pdbs false \
        -score:weights ref2015_cart \
        -fa_max_dis 9.0 \
        -ddg:legacy true \
        > "$log_file" 2>&1

    # Copy result
    if [ -f "mutations.ddg" ]; then
        cp "mutations.ddg" "$out_ddg"
        echo "[DONE] $basename ($(date +%H:%M:%S))"
    else
        echo "[FAIL] $basename - check $log_file"
    fi

    # Cleanup
    rm -rf "$tmpdir"
}

export -f run_one
export ROSETTA_BIN STRUCTURE OUT_DIR LOG_DIR

# Run in parallel
if command -v parallel &>/dev/null; then
    echo "Using GNU parallel ($N_JOBS jobs)..."
    printf '%s\n' "${MUT_FILES[@]}" | parallel -j "$N_JOBS" run_one {}
else
    echo "GNU parallel not found, using xargs ($N_JOBS jobs)..."
    printf '%s\n' "${MUT_FILES[@]}" | xargs -P "$N_JOBS" -I {} bash -c 'run_one "$@"' _ {}
fi

# Summary
echo ""
echo "========================================"
DONE_COUNT=$(ls "$OUT_DIR"/*.ddg 2>/dev/null | wc -l)
echo "Completed: $DONE_COUNT / $N_FILES"
echo "Results in: $OUT_DIR/"
echo "Logs in: $LOG_DIR/"
echo "========================================"
echo ""
echo "Next step: python 12_analyze_saturation.py"
