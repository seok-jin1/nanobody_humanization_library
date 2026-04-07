#!/usr/bin/env python3
"""
WT 및 돌연변이 MM-PBSA 안정성 결과 요약 (단순 버전)

Summarizes gmx_MMPBSA stability results, calculates DDG = DG(mut) - DG(WT).

Input:
    - mmpbsa_WT/FINAL_RESULTS.dat : WT MMPBSA result
    - mmpbsa_<MUT>/FINAL_RESULTS.dat : 13 mutant MMPBSA results
Output:
    - 02_mmpbsa_summary.txt : CSV format summary (Mutation, Position, DDG)

Usage:
    python 02_summarize_mmpbsa.py
"""

import os
import re
import glob

def parse_mmpbsa_result(result_file):
    """Parse gmx_MMPBSA FINAL_RESULTS.dat for total energy."""
    if not os.path.exists(result_file):
        return None
    with open(result_file) as f:
        for line in f:
            if 'TOTAL' in line and 'Std' not in line:
                parts = line.split()
                try:
                    return float(parts[1])  # Mean total energy
                except (IndexError, ValueError):
                    continue
    return None

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    wt_energy = parse_mmpbsa_result(f"{base_dir}/mmpbsa_WT/FINAL_RESULTS.dat")
    if wt_energy is None:
        print("WT results not found!")
        return

    mutations = {
        "Q1E": 1, "Q5V": 5, "S11L": 11, "A14P": 14,
        "G35S": 35, "V40A": 40, "A49S": 49, "I50V": 50,
        "A74S": 74, "K86R": 86, "P87A": 87, "M92V": 92, "Q120L": 120
    }

    print("=" * 70)
    print("MMPBSA MUTATION STABILITY ANALYSIS (ΔΔG in kcal/mol)")
    print("=" * 70)
    print(f"WT total energy: {wt_energy:.2f} kcal/mol")
    print()
    print(f"{'Mutation':<10} {'Pos':>4} {'ΔΔG (kcal/mol)':>16} {'Classification'}")
    print("-" * 60)

    results = []
    for mut, pos in sorted(mutations.items(), key=lambda x: x[1]):
        mut_energy = parse_mmpbsa_result(f"{base_dir}/mmpbsa_{mut}/FINAL_RESULTS.dat")
        if mut_energy is not None:
            ddg = mut_energy - wt_energy
            results.append((mut, pos, ddg))

    results_sorted = sorted(results, key=lambda x: x[2], reverse=True)

    for mut, pos, ddg in results_sorted:
        if ddg > 2.0:
            label = "❌ Strongly destabilizing"
        elif ddg > 1.0:
            label = "⚠️  Moderately destabilizing"
        elif ddg > 0.5:
            label = "⚡ Mildly destabilizing"
        elif ddg > -0.5:
            label = "✓  Neutral"
        else:
            label = "✅ Stabilizing"
        print(f"{mut:<10} {pos:>4} {ddg:>+16.2f}    {label}")

    print()
    print("=" * 70)

    # Save to file
    with open(f"{base_dir}/02_mmpbsa_summary.txt", "w") as f:
        f.write("Mutation,Position,DDG_kcal_mol\n")
        for mut, pos, ddg in results:
            f.write(f"{mut},{pos},{ddg:.3f}\n")
    print("Results saved to: 02_mmpbsa_summary.txt")

if __name__ == "__main__":
    main()
