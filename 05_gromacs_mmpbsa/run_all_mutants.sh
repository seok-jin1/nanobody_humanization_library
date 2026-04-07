#!/bin/bash
# Orchestrate full pipeline for all 13 mutants (sequential)
# Run AFTER: (1) 100 ns WT MD complete, (2) WT MMPBSA complete
#
# Usage: bash run_all_mutants.sh

set -e

BASE_DIR=/home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa
LOG="${BASE_DIR}/all_mutants.log"

MUTATIONS=(
    Q1E Q5V S11L A14P G35S V40A
    A49S I50V A74S K86R P87A M92V Q120L
)

echo "[$(date)] === Starting all 13 mutant pipelines ===" | tee "$LOG"
echo "Total mutations: ${#MUTATIONS[@]}" | tee -a "$LOG"

# ── Step 0: Generate mutant PDB files ────────────────────────────
echo "[$(date)] Generating mutant PDB files..." | tee -a "$LOG"
conda run -n docking-md python3 "${BASE_DIR}/make_mutants.py" 2>&1 | tee -a "$LOG"

# ── Step 1: Run pipeline for each mutant ─────────────────────────
TOTAL=${#MUTATIONS[@]}
COUNT=0
FAILED=()

for MUT in "${MUTATIONS[@]}"; do
    COUNT=$((COUNT + 1))
    echo "" | tee -a "$LOG"
    echo "[$(date)] [$COUNT/$TOTAL] === $MUT ===" | tee -a "$LOG"

    if bash "${BASE_DIR}/run_single_mutant.sh" "$MUT" 2>&1 | tee -a "$LOG"; then
        echo "[$(date)] [$COUNT/$TOTAL] $MUT: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] [$COUNT/$TOTAL] $MUT: FAILED" | tee -a "$LOG"
        FAILED+=("$MUT")
    fi
done

# ── Step 2: Summarize results ─────────────────────────────────────
echo "" | tee -a "$LOG"
echo "[$(date)] === All pipelines complete ===" | tee -a "$LOG"
echo "Failed: ${FAILED[*]:-none}" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "[$(date)] Generating summary..." | tee -a "$LOG"
python3 "${BASE_DIR}/summarize_all_mmpbsa.py" 2>&1 | tee -a "$LOG"
