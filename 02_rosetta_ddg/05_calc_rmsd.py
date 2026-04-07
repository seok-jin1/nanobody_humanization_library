#!/usr/bin/env python3
"""
AlphaFold3 모델과 템플릿 모델 간 Kabsch 정렬 RMSD 계산

Calculates the C-alpha RMSD between an AlphaFold3-predicted nanobody
structure (mmCIF) and the 1ZVH template-based model (PDB) using the
Kabsch algorithm for optimal superposition. Reports per-residue
deviations to identify regions of structural divergence.

Input:
    - ../01_structure_prediction/fold_anti_fap_nb_model_0.cif : AF3 model (mmCIF)
    - relaxed/threaded_template_1ZVH_chainL_0001_0005.pdb : template model (PDB)
Output:
    - console : RMSD value, per-residue deviations, reliability assessment

Usage:
    python 05_calc_rmsd.py
"""

import numpy as np
from collections import defaultdict

def parse_pdb_coords(pdb_file):
    """Extract CA coordinates from PDB file."""
    coords = []
    residues = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resnum = int(line[22:26].strip())
                coords.append([x, y, z])
                residues.append(resnum)

    return np.array(coords), residues

def parse_cif_coords(cif_file):
    """Extract CA coordinates from mmCIF file."""
    coords = []
    residues = []

    with open(cif_file, 'r') as f:
        for line in f:
            # Look for ATOM lines with CA
            if line.startswith('ATOM') and ' CA ' in line:
                parts = line.split()
                if len(parts) < 18:
                    continue
                try:
                    # Format: ATOM id type label_atom_id . comp asym entity seq ? x y z occ B auth_seq auth_asym model
                    # Example: ATOM 2   C CA  . GLN A 1 1   ? 8.197   -0.517  19.493  1.00 91.27 1   A 1
                    x = float(parts[10])  # Cartn_x
                    y = float(parts[11])  # Cartn_y
                    z = float(parts[12])  # Cartn_z
                    resnum = int(parts[8])  # label_seq_id
                    coords.append([x, y, z])
                    residues.append(resnum)
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line: {line.strip()}")
                    print(f"Error: {e}")
                    continue

    return np.array(coords), residues

def kabsch_rmsd(coords1, coords2):
    """Calculate RMSD after optimal rotation using Kabsch algorithm."""
    if len(coords1) != len(coords2):
        min_len = min(len(coords1), len(coords2))
        print(f"Warning: Different lengths ({len(coords1)} vs {len(coords2)}), using first {min_len}")
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]

    # Center coordinates
    center1 = np.mean(coords1, axis=0)
    center2 = np.mean(coords2, axis=0)
    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2

    # Kabsch algorithm - find optimal rotation
    # Compute covariance matrix
    H = coords1_centered.T @ coords2_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute optimal rotation
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T

    # Apply rotation
    coords1_rotated = coords1_centered @ R

    # Calculate RMSD
    diff = coords1_rotated - coords2_centered
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd, coords1_rotated + center2

def main():
    print("="*80)
    print("RMSD CALCULATION: AlphaFold3 vs Template-Based Model")
    print("="*80)
    print()

    af3_file = "../01_structure_prediction/fold_anti_fap_nb_model_0.cif"
    template_file = "relaxed/threaded_template_1ZVH_chainL_0001_0005.pdb"

    print(f"AlphaFold3 model: {af3_file}")
    print(f"Template model:   {template_file}")
    print()

    # Parse coordinates
    print("Parsing structures...")
    af3_coords, af3_res = parse_cif_coords(af3_file)
    template_coords, template_res = parse_pdb_coords(template_file)

    print(f"  AF3:      {len(af3_coords)} CA atoms")
    print(f"  Template: {len(template_coords)} CA atoms")
    print()

    if len(af3_coords) == 0 or len(template_coords) == 0:
        print("ERROR: Could not parse coordinates. Check file formats.")
        return

    # Calculate RMSD with Kabsch alignment
    print("Calculating RMSD with optimal superposition (Kabsch)...")
    rmsd, af3_aligned = kabsch_rmsd(af3_coords, template_coords)

    print()
    print("="*80)
    print("RESULTS")
    print("="*80)
    print(f"C-alpha RMSD (Kabsch-aligned): {rmsd:.2f} Å")
    print()

    if rmsd < 2.0:
        print("✅ EXCELLENT: Structures are very similar (RMSD < 2 Å)")
        print("   → ddG calculations on template model are reliable")
    elif rmsd < 3.0:
        print("✓ GOOD: Structures are similar (RMSD < 3 Å)")
        print("   → ddG calculations should be consistent")
    elif rmsd < 5.0:
        print("⚠ MODERATE: Some structural differences (RMSD 3-5 Å)")
        print("   → Consider comparing ddG results from both structures")
    else:
        print("❌ HIGH: Significant structural differences (RMSD > 5 Å)")
        print("   → Recommend re-running ddG with AlphaFold3 structure")

    print()
    print("NOTE: This is a Kabsch-aligned RMSD (optimal superposition).")
    print()

    # Calculate per-residue distance (using aligned structures)
    min_len = min(len(af3_aligned), len(template_coords))
    if min_len > 0:
        print("="*80)
        print("PER-RESIDUE DEVIATIONS (largest differences)")
        print("="*80)

        deviations = np.sqrt(np.sum((af3_aligned[:min_len] - template_coords[:min_len])**2, axis=1))
        sorted_indices = np.argsort(deviations)[::-1]

        print("\nTop 10 most different residues:")
        print(f"{'Residue':>8} {'Deviation (Å)':>15}")
        print("-"*25)
        for i in sorted_indices[:10]:
            print(f"{af3_res[i]:8d} {deviations[i]:15.2f}")

        print("\nThese differences may indicate:")
        print("  - CDR loop conformational differences")
        print("  - Side-chain packing differences")
        print("  - Local structural variations")

    print()
    print("="*80)

if __name__ == "__main__":
    main()
