#!/bin/bash
set -e
GMX=/usr/local/gromacs/bin/gmx   # system GROMACS with CUDA support

cd /home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa

echo "[$(date)] === NVT equilibration (100 ps) ==="
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu 2>&1
echo "[$(date)] NVT done"

echo "[$(date)] === NPT equilibration (100 ps) ==="
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu 2>&1
echo "[$(date)] NPT done"

echo "[$(date)] === Production MD (100 ns) ==="
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2 2>&1
$GMX mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu 2>&1
echo "[$(date)] Production MD done"

echo "[$(date)] === All MD complete! ==="
