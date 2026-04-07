#!/usr/bin/env python3
"""
WT 및 돌연변이 MM-PBSA 안정성 결과 요약 (오차 전파 포함)

Enhanced MMPBSA summary with error propagation (std dev tracking).
DDG = DG(mutant) - DG(WT), with propagated uncertainty.

Input:
    - mmpbsa_WT/FINAL_RESULTS.dat : WT MMPBSA result (last 10 ns of 100 ns trajectory)
    - mutants/<MUT>/mmpbsa/FINAL_RESULTS.dat : 13 mutant MMPBSA results (10 ns each)
Output:
    - 03_mmpbsa_ddg_summary.csv : CSV with DeltaG, Std, DeltaDeltaG, Error

Usage:
    python 03_summarize_all_mmpbsa.py
"""

import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

MUTATIONS = [
    ("Q1E", 1), ("Q5V", 5), ("S11L", 11), ("A14P", 14),
    ("G35S", 35), ("V40A", 40), ("A49S", 49), ("I50V", 50),
    ("A74S", 74), ("K86R", 86), ("P87A", 87), ("M92V", 92), ("Q120L", 120),
]


def parse_final_results(dat_file):
    """Parse FINAL_RESULTS.dat → return (mean_total, std_total) or None."""
    if not os.path.exists(dat_file):
        return None
    mean = std = None
    with open(dat_file) as f:
        for line in f:
            if line.strip().startswith("TOTAL") and "Std" not in line:
                parts = line.split()
                try:
                    mean = float(parts[1])
                    std = float(parts[2]) if len(parts) > 2 else 0.0
                    return mean, std
                except (IndexError, ValueError):
                    continue
    return None


def main():
    # ── WT ──────────────────────────────────────────────────────
    wt_file = os.path.join(BASE_DIR, "mmpbsa_WT", "FINAL_RESULTS.dat")
    wt_result = parse_final_results(wt_file)
    if wt_result is None:
        print("ERROR: WT results not found at", wt_file)
        return
    wt_mean, wt_std = wt_result

    # ── Mutants ─────────────────────────────────────────────────
    results = []
    missing = []
    for mut_name, pos in MUTATIONS:
        dat_file = os.path.join(BASE_DIR, "mutants", mut_name, "mmpbsa", "FINAL_RESULTS.dat")
        result = parse_final_results(dat_file)
        if result is None:
            missing.append(mut_name)
            continue
        mut_mean, mut_std = result
        ddg = mut_mean - wt_mean
        ddg_err = (wt_std**2 + mut_std**2) ** 0.5
        results.append((mut_name, pos, ddg, ddg_err, mut_mean, mut_std))

    # ── Print ────────────────────────────────────────────────────
    print("=" * 72)
    print("MMPBSA MUTATION STABILITY ANALYSIS (ΔΔG in kcal/mol)")
    print("WT: last 10 ns of 100 ns MD | Mutants: 10 ns MD each")
    print("=" * 72)
    print(f"WT ΔG(stability): {wt_mean:.2f} ± {wt_std:.2f} kcal/mol")
    print()
    print(f"{'Mutation':<10} {'Pos':>4} {'ΔΔG':>10} {'±':>8}  {'Classification'}")
    print("-" * 72)

    for mut_name, pos, ddg, ddg_err, mut_mean, mut_std in sorted(results, key=lambda x: x[2], reverse=True):
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
        print(f"{mut_name:<10} {pos:>4} {ddg:>+10.2f} {ddg_err:>8.2f}  {label}")

    if missing:
        print()
        print(f"Missing results: {', '.join(missing)}")

    print()
    print("=" * 72)

    # ── Save CSV ─────────────────────────────────────────────────
    out_csv = os.path.join(BASE_DIR, "03_mmpbsa_ddg_summary.csv")
    with open(out_csv, "w") as f:
        f.write("Mutation,Position,DeltaG_mut,Std_mut,DeltaDeltaG,Error\n")
        for mut_name, pos, ddg, ddg_err, mut_mean, mut_std in results:
            f.write(f"{mut_name},{pos},{mut_mean:.3f},{mut_std:.3f},{ddg:.3f},{ddg_err:.3f}\n")
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
