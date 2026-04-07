#!/usr/bin/env python3
"""
WT vs Humanized MD trajectory 비교 분석

GROMACS MD 결과에서 RMSD, RMSF, Rg를 추출하여
WT와 Humanized nanobody의 구조적 안정성을 비교한다.

Input:
    - WT MD: ../05_gromacs_mmpbsa/ (또는 지정 경로)의 xvg 파일
    - Humanized MD: md_humanized/의 xvg 파일
    이 스크립트 실행 전에 GROMACS로 xvg 파일을 먼저 생성해야 함 (03_extract_xvg.sh)
Output:
    - 03_md_comparison.csv                   : 비교 통계 데이터
    - 03_md_comparison.txt                   : 텍스트 보고서
    - figures/03_rmsd_comparison.png/pdf     : RMSD 시계열 비교
    - figures/03_rmsf_comparison.png/pdf     : RMSF per-residue 비교
    - figures/03_rg_comparison.png/pdf       : Radius of gyration 비교

Usage:
    # 1. xvg 파일 추출
    bash 03_extract_xvg.sh
    # 2. 비교 분석
    python 03_compare_md.py
"""

import csv
import os

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

# Paths to xvg files — adjust if WT trajectory is in a different location
WT_DIR = os.path.join(SCRIPT_DIR, "md_wt_xvg")
HUM_DIR = os.path.join(SCRIPT_DIR, "md_humanized")


def parse_xvg(filepath):
    """Parse GROMACS .xvg file → (x_array, y_array)."""
    x, y = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)


def plot_timeseries(wt_file, hum_file, title, ylabel, out_name):
    """Plot WT vs Humanized time series comparison."""
    if not os.path.exists(wt_file) or not os.path.exists(hum_file):
        print(f"  Skipping {out_name}: file not found")
        return None, None

    wt_x, wt_y = parse_xvg(wt_file)
    hum_x, hum_y = parse_xvg(hum_file)

    # Convert ps → ns if needed
    if wt_x[-1] > 10000:
        wt_x = wt_x / 1000
        hum_x = hum_x / 1000

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(wt_x, wt_y, color='#ef5350', alpha=0.7, linewidth=0.5, label='WT')
    ax.plot(hum_x, hum_y, color='#42a5f5', alpha=0.7, linewidth=0.5, label='Humanized')

    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, f"{out_name}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGURES_DIR, f"{out_name}.pdf"), bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return (np.mean(wt_y), np.std(wt_y)), (np.mean(hum_y), np.std(hum_y))


def plot_rmsf(wt_file, hum_file):
    """Plot per-residue RMSF comparison."""
    if not os.path.exists(wt_file) or not os.path.exists(hum_file):
        print("  Skipping RMSF: file not found")
        return None, None

    wt_x, wt_y = parse_xvg(wt_file)
    hum_x, hum_y = parse_xvg(hum_file)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(wt_x, wt_y, color='#ef5350', alpha=0.8, linewidth=1, label='WT')
    ax.plot(hum_x, hum_y, color='#42a5f5', alpha=0.8, linewidth=1, label='Humanized')

    # Highlight mutation positions (PDB numbering)
    mut_positions = [5, 11, 35, 50, 74, 86, 92, 120]
    for pos in mut_positions:
        ax.axvline(pos, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('Residue Number', fontsize=11)
    ax.set_ylabel('RMSF (nm)', fontsize=11)
    ax.set_title('Per-Residue RMSF: WT vs Humanized', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, "03_rmsf_comparison.png"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(FIGURES_DIR, "03_rmsf_comparison.pdf"), bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return (np.mean(wt_y), np.std(wt_y)), (np.mean(hum_y), np.std(hum_y))


def main():
    print("WT vs Humanized MD Comparison")
    print("=" * 50)

    stats = {}

    # RMSD
    rmsd = plot_timeseries(
        os.path.join(WT_DIR, "rmsd.xvg"),
        os.path.join(HUM_DIR, "rmsd.xvg"),
        "Backbone RMSD: WT vs Humanized",
        "RMSD (nm)",
        "03_rmsd_comparison"
    )
    if rmsd[0]:
        stats['rmsd'] = rmsd
        print(f"  RMSD — WT: {rmsd[0][0]:.3f}±{rmsd[0][1]:.3f} nm, Hum: {rmsd[1][0]:.3f}±{rmsd[1][1]:.3f} nm")

    # Rg
    rg = plot_timeseries(
        os.path.join(WT_DIR, "gyrate.xvg"),
        os.path.join(HUM_DIR, "gyrate.xvg"),
        "Radius of Gyration: WT vs Humanized",
        "Rg (nm)",
        "03_rg_comparison"
    )
    if rg[0]:
        stats['rg'] = rg
        print(f"  Rg   — WT: {rg[0][0]:.3f}±{rg[0][1]:.3f} nm, Hum: {rg[1][0]:.3f}±{rg[1][1]:.3f} nm")

    # RMSF
    rmsf = plot_rmsf(
        os.path.join(WT_DIR, "rmsf.xvg"),
        os.path.join(HUM_DIR, "rmsf.xvg"),
    )
    if rmsf[0]:
        stats['rmsf'] = rmsf
        print(f"  RMSF — WT: {rmsf[0][0]:.3f}±{rmsf[0][1]:.3f} nm, Hum: {rmsf[1][0]:.3f}±{rmsf[1][1]:.3f} nm")

    if not stats:
        print("\nNo xvg files found. Run 03_extract_xvg.sh first.")
        return

    # CSV
    out_csv = os.path.join(SCRIPT_DIR, "03_md_comparison.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "wt_mean", "wt_std", "hum_mean", "hum_std", "unit"])
        for metric, (wt_stat, hum_stat) in stats.items():
            unit = "nm"
            writer.writerow([metric, f"{wt_stat[0]:.4f}", f"{wt_stat[1]:.4f}",
                           f"{hum_stat[0]:.4f}", f"{hum_stat[1]:.4f}", unit])
    print(f"\nSaved: {out_csv}")

    # Text report
    out_txt = os.path.join(SCRIPT_DIR, "03_md_comparison.txt")
    with open(out_txt, "w") as f:
        f.write("WT vs Humanized MD Comparison\n")
        f.write("=" * 50 + "\n\n")
        for metric, (wt_stat, hum_stat) in stats.items():
            f.write(f"{metric.upper()}:\n")
            f.write(f"  WT:        {wt_stat[0]:.4f} ± {wt_stat[1]:.4f} nm\n")
            f.write(f"  Humanized: {hum_stat[0]:.4f} ± {hum_stat[1]:.4f} nm\n\n")
    print(f"Saved: {out_txt}")


if __name__ == "__main__":
    main()
