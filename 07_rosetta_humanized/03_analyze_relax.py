#!/usr/bin/env python3
"""
WT vs Humanized FastRelax 결과 비교 분석 및 시각화

FastRelax로 생성된 score 파일을 파싱하여 total_score를 비교하고,
multi-mutant ddG 결과도 함께 보고한다.

Input:
    - 01_multi_mutant.ddg     : multi-mutant ddG raw output
    - 02_relax_wt.sc          : WT relaxed scores (5 structures)
    - 02_relax_humanized.sc   : Humanized relaxed scores (5 structures)
Output:
    - 03_stability_comparison.csv          : 비교 데이터
    - 03_stability_comparison.txt          : 텍스트 보고서
    - figures/03_stability_comparison.png  : total score 비교 box plot
    - figures/03_stability_comparison.pdf  : PDF

Usage:
    python 03_analyze_relax.py
"""

import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

DDG_FILE = os.path.join(SCRIPT_DIR, "01_multi_mutant.ddg")
WT_SC = os.path.join(SCRIPT_DIR, "02_relax_wt.sc")
HUM_SC = os.path.join(SCRIPT_DIR, "02_relax_humanized.sc")


def parse_scorefile(filepath):
    """Parse Rosetta score file → list of total_score values."""
    scores = []
    if not os.path.exists(filepath):
        return scores
    with open(filepath) as f:
        header = None
        for line in f:
            if line.startswith("SEQUENCE:"):
                continue
            if line.startswith("SCORE:") and header is None:
                header = line.split()
                continue
            if line.startswith("SCORE:") and header is not None:
                parts = line.split()
                try:
                    idx = header.index("total_score")
                    scores.append(float(parts[idx]))
                except (ValueError, IndexError):
                    continue
    return scores


def parse_ddg_file(filepath):
    """Parse multi-mutant .ddg file → (mean_ddg, std_ddg)."""
    if not os.path.exists(filepath):
        return None, None
    wt_energies = []
    mut_energies = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith("COMPLEX:"):
                continue
            parts = line.split()
            name = parts[2].rstrip(':')
            energy = float(parts[3])
            if name == 'WT':
                wt_energies.append(energy)
            else:
                mut_energies.append(energy)

    if not wt_energies or not mut_energies:
        return None, None
    ddg = np.mean(mut_energies) - np.mean(wt_energies)
    std = np.sqrt(np.std(mut_energies, ddof=1)**2 + np.std(wt_energies, ddof=1)**2) if len(mut_energies) > 1 else 0
    return ddg, std


def main():
    wt_scores = parse_scorefile(WT_SC)
    hum_scores = parse_scorefile(HUM_SC)
    ddg, ddg_std = parse_ddg_file(DDG_FILE)

    has_relax = bool(wt_scores and hum_scores)
    has_ddg = ddg is not None

    if not has_relax and not has_ddg:
        print("No results found. Run 01_run_multi_mutant_ddg.sh and/or 02_run_relax_compare.sh first.")
        return

    # --- Report ---
    lines = []
    lines.append("WT vs Humanized Stability Comparison")
    lines.append("=" * 50)
    lines.append("")

    if has_ddg:
        lines.append(f"1. Multi-Mutant Cartesian ddG (8 simultaneous mutations)")
        lines.append(f"   ddG = {ddg:+.3f} ± {ddg_std:.3f} REU")
        if ddg < -0.5:
            lines.append(f"   → Stabilizing (humanization improves stability)")
        elif ddg > 2.0:
            lines.append(f"   → Destabilizing (humanization harms stability)")
        else:
            lines.append(f"   → Neutral (minimal impact)")
        lines.append("")

    if has_relax:
        wt_mean = np.mean(wt_scores)
        wt_std = np.std(wt_scores, ddof=1)
        hum_mean = np.mean(hum_scores)
        hum_std = np.std(hum_scores, ddof=1)
        diff = hum_mean - wt_mean

        lines.append(f"2. FastRelax Total Score Comparison")
        lines.append(f"   WT:        {wt_mean:.2f} ± {wt_std:.2f} REU (n={len(wt_scores)})")
        lines.append(f"   Humanized: {hum_mean:.2f} ± {hum_std:.2f} REU (n={len(hum_scores)})")
        lines.append(f"   Difference: {diff:+.2f} REU")
        lines.append(f"   → {'Humanized more stable' if diff < 0 else 'WT more stable' if diff > 0 else 'Similar stability'}")
        lines.append("")

    report = '\n'.join(lines)
    print(report)

    out_txt = os.path.join(SCRIPT_DIR, "03_stability_comparison.txt")
    with open(out_txt, "w") as f:
        f.write(report)
    print(f"Saved: {out_txt}")

    # --- CSV ---
    out_csv = os.path.join(SCRIPT_DIR, "03_stability_comparison.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "wt_value", "humanized_value", "difference", "unit"])
        if has_ddg:
            writer.writerow(["multi_mutant_ddg", "0", f"{ddg:.3f}", f"{ddg:.3f}", "REU"])
        if has_relax:
            writer.writerow(["relax_total_score", f"{wt_mean:.2f}", f"{hum_mean:.2f}",
                            f"{diff:.2f}", "REU"])
    print(f"Saved: {out_csv}")

    # --- Visualization ---
    if has_relax:
        fig, ax = plt.subplots(figsize=(5, 5))

        bp = ax.boxplot([wt_scores, hum_scores], labels=['WT', 'Humanized'],
                        patch_artist=True, widths=0.5)
        bp['boxes'][0].set_facecolor('#ef5350')
        bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][1].set_facecolor('#42a5f5')
        bp['boxes'][1].set_alpha(0.7)

        # Individual points
        for i, scores in enumerate([wt_scores, hum_scores]):
            jitter = np.random.normal(0, 0.04, len(scores))
            ax.scatter([i + 1 + j for j in jitter], scores, color='black',
                      s=20, alpha=0.6, zorder=3)

        ax.set_ylabel('Total Score (REU)', fontsize=12, fontweight='bold')
        ax.set_title('FastRelax Total Score: WT vs Humanized',
                    fontsize=13, fontweight='bold', pad=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, linestyle='--', alpha=0.2)

        if has_ddg:
            ax.text(0.02, 0.98, f"Multi-mutant ΔΔG = {ddg:+.2f} ± {ddg_std:.2f} REU",
                   transform=ax.transAxes, fontsize=9, va='top',
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                            edgecolor='gray', alpha=0.8))

        fig.tight_layout()

        out_png = os.path.join(FIGURES_DIR, "03_stability_comparison.png")
        out_pdf = os.path.join(FIGURES_DIR, "03_stability_comparison.pdf")
        fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
        print(f"Saved: {out_png}")
        print(f"Saved: {out_pdf}")
        plt.close(fig)


if __name__ == "__main__":
    main()
