#!/usr/bin/env python3
"""
13개 humanization 돌연변이체 PDB 파일 생성

Generates 13 mutant PDB files from WT nanobody using pdbfixer/OpenMM.
Each mutant gets its own directory under mutants/.

Input:
    - ../01_structure_prediction/nanobody_af3.pdb : WT nanobody structure
Output:
    - mutants/<MUT_NAME>/nanobody_<MUT_NAME>.pdb : 13 mutant PDB files

Usage:
    python 01_make_mutants.py

Dependencies:
    pdbfixer, OpenMM
"""

import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WT_PDB = os.path.join(BASE_DIR, "..", "01_structure_prediction", "nanobody_af3.pdb")
MUTANTS_DIR = os.path.join(BASE_DIR, "mutants")

# 13 humanization mutations
# Format: "name": ("pdbfixer_string", sequential_resnum)
# pdbfixer format: "ORIGRES-RESNUM-NEWRES" (3-letter codes)
MUTATIONS = {
    "Q1E":   "GLN-1-GLU",
    "Q5V":   "GLN-5-VAL",
    "S11L":  "SER-11-LEU",
    "A14P":  "ALA-14-PRO",
    "G35S":  "GLY-35-SER",
    "V40A":  "VAL-40-ALA",
    "A49S":  "ALA-49-SER",
    "I50V":  "ILE-50-VAL",
    "A74S":  "ALA-74-SER",
    "K86R":  "LYS-86-ARG",
    "P87A":  "PRO-87-ALA",
    "M92V":  "MET-92-VAL",
    "Q120L": "GLN-120-LEU",
}

os.makedirs(MUTANTS_DIR, exist_ok=True)

success = []
failed = []

for mut_name, mut_string in MUTATIONS.items():
    out_dir = os.path.join(MUTANTS_DIR, mut_name)
    os.makedirs(out_dir, exist_ok=True)
    out_pdb = os.path.join(out_dir, f"nanobody_{mut_name}.pdb")

    print(f"[{mut_name}] Applying mutation: {mut_string} ...", end=" ", flush=True)
    try:
        fixer = PDBFixer(filename=WT_PDB)
        fixer.applyMutations([mut_string], "A")
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        # Do NOT add hydrogens here — pdb2gmx will handle them with -ignh

        with open(out_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        # Verify atom count
        with open(out_pdb) as f:
            atoms = sum(1 for l in f if l.startswith("ATOM"))
        print(f"OK ({atoms} heavy atoms) → {out_pdb}")
        success.append(mut_name)

    except Exception as e:
        print(f"FAILED: {e}")
        failed.append((mut_name, str(e)))

print()
print("=" * 50)
print(f"Success: {len(success)}/{len(MUTATIONS)}")
if failed:
    print("Failed:")
    for name, err in failed:
        print(f"  {name}: {err}")
else:
    print("All mutations generated successfully!")
print(f"\nOutput directory: {MUTANTS_DIR}/")
