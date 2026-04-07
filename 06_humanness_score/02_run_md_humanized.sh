#!/bin/bash
# Humanized nanobody MD simulation (GROMACS)
#
# WT와 동일한 조건으로 humanized nanobody의 100ns MD를 수행.
# 05_gromacs_mmpbsa/의 mdp 파일을 공유하며, 별도 작업 디렉토리에서 실행.
#
# Input:
#   - ../01_structure_prediction/fold_humanized_fap_nb_model_0.cif : Humanized AF3 structure
#   - ../05_gromacs_mmpbsa/*.mdp : MD parameter files
# Output:
#   - md_humanized/ : GROMACS trajectory + analysis files
#
# Prerequisites:
#   - GROMACS with GPU support
#   - gemmi (CIF → PDB conversion): pip install gemmi
#
# Usage:
#   bash 02_run_md_humanized.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GMX=/usr/local/gromacs/bin/gmx
WORK_DIR="${SCRIPT_DIR}/md_humanized"
MDP_DIR="${SCRIPT_DIR}/../05_gromacs_mmpbsa"
STRUCTURE="${SCRIPT_DIR}/../01_structure_prediction/fold_humanized_fap_nb_model_0.cif"

mkdir -p "$WORK_DIR"

if [ ! -f "$STRUCTURE" ]; then
    echo "ERROR: Humanized structure not found: $STRUCTURE"
    exit 1
fi

echo "========================================"
echo "Humanized Nanobody MD Simulation"
echo "========================================"
echo "Structure: $STRUCTURE"
echo "Work dir:  $WORK_DIR"
echo ""

# Convert CIF → PDB
echo "[$(date +%H:%M:%S)] Converting CIF → PDB..."
python3 -c "
import gemmi
st = gemmi.read_structure('$STRUCTURE')
st.write_pdb('$WORK_DIR/humanized.pdb')
print('Converted to PDB')
"

cd "$WORK_DIR"

# Generate topology
echo "[$(date +%H:%M:%S)] Generating topology (AMBER99SB-ILDN)..."
echo "1" | $GMX pdb2gmx -f humanized.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh 2>&1

# Solvation
echo "[$(date +%H:%M:%S)] Solvation..."
$GMX editconf -f processed.gro -o box.gro -c -d 1.2 -bt dodecahedron 2>&1
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top 2>&1

# Ions
echo "[$(date +%H:%M:%S)] Adding ions..."
$GMX grompp -f "$MDP_DIR/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1
echo "SOL" | $GMX genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral 2>&1

# Energy minimization
echo "[$(date +%H:%M:%S)] Energy minimization..."
$GMX grompp -f "$MDP_DIR/em.mdp" -c ions.gro -p topol.top -o em.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm em -ntmpi 1 -ntomp 8 2>&1
echo "[$(date +%H:%M:%S)] EM done"

# NVT equilibration
echo "[$(date +%H:%M:%S)] NVT equilibration (100 ps)..."
$GMX grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu 2>&1
echo "[$(date +%H:%M:%S)] NVT done"

# NPT equilibration
echo "[$(date +%H:%M:%S)] NPT equilibration (100 ps)..."
$GMX grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu 2>&1
echo "[$(date +%H:%M:%S)] NPT done"

# Production MD (100 ns)
echo "[$(date +%H:%M:%S)] Production MD (100 ns)..."
$GMX grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu 2>&1
echo "[$(date +%H:%M:%S)] Production MD done"

echo ""
echo "========================================"
echo "MD complete. Trajectory: $WORK_DIR/md.xtc"
echo "Next: python 03_compare_md.py"
echo "========================================"
