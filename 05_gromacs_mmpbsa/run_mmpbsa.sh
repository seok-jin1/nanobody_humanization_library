#!/bin/bash
# gmx_MMPBSA WT Stability Analysis
# Run after MD production is complete

set -e
GMX=/usr/local/gromacs/bin/gmx   # system GROMACS 2024.1 (md.tpr 생성 버전)

source /home/laugh/miniconda3/etc/profile.d/conda.sh
conda activate docking-md

cd /home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa

echo "[$(date)] Preparing trajectory for MMPBSA..."

# Center protein and remove PBC artifacts
$GMX trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center 2>&1 <<< $'1\n0'

echo "[$(date)] Running WT stability calculation..."
mkdir -p mmpbsa_WT

gmx_MMPBSA \
    -O \
    -s \
    -i mmpbsa_stability.in \
    -cs md.tpr \
    -ct md_noPBC.xtc \
    -cp topol.top \
    -o mmpbsa_WT/FINAL_RESULTS.dat \
    -eo mmpbsa_WT/FINAL_RESULTS.csv \
    --prefix mmpbsa_WT/ \
    -nogui \
    2>&1 | tee mmpbsa_WT/run.log

echo "[$(date)] === WT MMPBSA complete! ==="
echo "Results: mmpbsa_WT/FINAL_RESULTS.dat"
