#!/bin/bash
# Full pipeline for a single mutant: EM → NVT → NPT → 10 ns MD → MMPBSA
# Usage: bash run_single_mutant.sh <MUT_NAME>
# Example: bash run_single_mutant.sh Q1E

set -e

MUT_NAME="$1"
if [ -z "$MUT_NAME" ]; then
    echo "Usage: $0 <MUT_NAME>"
    echo "Example: $0 Q1E"
    exit 1
fi

GMX=/usr/local/gromacs/bin/gmx
BASE_DIR=/home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa
MUT_DIR="${BASE_DIR}/mutants/${MUT_NAME}"
INPUT_PDB="${MUT_DIR}/nanobody_${MUT_NAME}.pdb"
LOG="${MUT_DIR}/pipeline.log"

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: $INPUT_PDB not found. Run make_mutants.py first."
    exit 1
fi

cd "$MUT_DIR"
echo "[$(date)] === Starting pipeline for $MUT_NAME ===" | tee -a "$LOG"

# ── Step 1: pdb2gmx ──────────────────────────────────────────────
echo "[$(date)] pdb2gmx..." | tee -a "$LOG"
$GMX pdb2gmx \
    -f "nanobody_${MUT_NAME}.pdb" \
    -o processed.gro \
    -p topol.top \
    -i posre.itp \
    -ignh \
    -ff amber99sb-ildn \
    -water tip3p \
    -maxwarn 2 >> "$LOG" 2>&1

# Fix posre.itp include path
sed -i 's|#include ".*posre.itp"|#include "posre.itp"|' topol.top

# ── Step 2: editconf ─────────────────────────────────────────────
echo "[$(date)] editconf..." | tee -a "$LOG"
$GMX editconf \
    -f processed.gro \
    -o box.gro \
    -c -d 1.2 -bt dodecahedron >> "$LOG" 2>&1

# ── Step 3: solvate ──────────────────────────────────────────────
echo "[$(date)] solvate..." | tee -a "$LOG"
$GMX solvate \
    -cp box.gro \
    -cs spc216.gro \
    -o solv.gro \
    -p topol.top >> "$LOG" 2>&1

# ── Step 4: genion ───────────────────────────────────────────────
echo "[$(date)] genion..." | tee -a "$LOG"
$GMX grompp \
    -f "${BASE_DIR}/ions.mdp" \
    -c solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 2 >> "$LOG" 2>&1
echo "SOL" | $GMX genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p topol.top \
    -pname NA -nname CL -neutral >> "$LOG" 2>&1

# ── Step 5: Energy minimization ──────────────────────────────────
echo "[$(date)] Energy minimization..." | tee -a "$LOG"
$GMX grompp \
    -f "${BASE_DIR}/em.mdp" \
    -c solv_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 2 >> "$LOG" 2>&1
$GMX mdrun -deffnm em -ntmpi 1 -ntomp 4 >> "$LOG" 2>&1

# ── Step 6: NVT equilibration (100 ps) ───────────────────────────
echo "[$(date)] NVT equilibration..." | tee -a "$LOG"
$GMX grompp \
    -f "${BASE_DIR}/nvt.mdp" \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr \
    -maxwarn 2 >> "$LOG" 2>&1
$GMX mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu >> "$LOG" 2>&1

# ── Step 7: NPT equilibration (100 ps) ───────────────────────────
echo "[$(date)] NPT equilibration..." | tee -a "$LOG"
$GMX grompp \
    -f "${BASE_DIR}/npt.mdp" \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt.tpr \
    -maxwarn 2 >> "$LOG" 2>&1
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu >> "$LOG" 2>&1

# ── Step 8: Production MD (10 ns) ────────────────────────────────
echo "[$(date)] Production MD (10 ns)..." | tee -a "$LOG"
$GMX grompp \
    -f "${BASE_DIR}/md_10ns.mdp" \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -o md.tpr \
    -maxwarn 2 >> "$LOG" 2>&1
$GMX mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu >> "$LOG" 2>&1

# ── Step 9: MMPBSA ───────────────────────────────────────────────
echo "[$(date)] gmx_MMPBSA..." | tee -a "$LOG"
source /home/laugh/miniconda3/etc/profile.d/conda.sh
conda activate docking-md

$GMX trjconv \
    -s md.tpr -f md.xtc \
    -o md_noPBC.xtc \
    -pbc mol -center >> "$LOG" 2>&1 <<< $'1\n0'

mkdir -p mmpbsa
gmx_MMPBSA \
    -O \
    -s \
    -i "${BASE_DIR}/mmpbsa_stability.in" \
    -cs md.tpr \
    -ct md_noPBC.xtc \
    -cp topol.top \
    -o mmpbsa/FINAL_RESULTS.dat \
    -eo mmpbsa/FINAL_RESULTS.csv \
    --prefix mmpbsa/ \
    -nogui >> "$LOG" 2>&1

echo "[$(date)] === $MUT_NAME pipeline complete! ===" | tee -a "$LOG"
echo "Results: ${MUT_DIR}/mmpbsa/FINAL_RESULTS.dat"
