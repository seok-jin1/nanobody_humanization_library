#!/usr/bin/env python3
"""
ddG 결과를 시각화하여 과제/논문용 그림 생성

Generates publication-ready bar plots of Cartesian ddG values:
  Panel A: Set1 - 13 humanization mutations sorted by IMGT position
  Panel B: Set2 - I55 and P96 saturation mutagenesis
Bars are colored by stability classification. IMGT positions and
FR/CDR regions are annotated.

Input:
    - 02_ddg_predictions_final.out : Set1 ddG predictions
    - 02_ddg_predictions_set2.out  : Set2 ddG predictions
Output:
    - figures/04_ddg_barplot.png : bar plot (PNG)
    - figures/04_ddg_barplot.pdf : bar plot (PDF)
    - 04_ddg_plot_data.csv       : underlying plot data

Usage:
    python 04_visualize_ddg.py
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

SET1_FILE = os.path.join(SCRIPT_DIR, "02_ddg_predictions_final.out")
SET2_FILE = os.path.join(SCRIPT_DIR, "02_ddg_predictions_set2.out")

# Rosetta MUT name → (IMGT pos, region, human-readable label)
MUTATION_INFO = {
    # Set1
    'MUT_1GLU':   (1,   'FR1',  'Q1E'),
    'MUT_5VAL':   (5,   'FR1',  'Q5V'),
    'MUT_11LEU':  (12,  'FR1',  'S12L'),
    'MUT_14PRO':  (15,  'FR1',  'A15P'),
    'MUT_35SER':  (40,  'FR2',  'G40S'),
    'MUT_40ALA':  (45,  'FR2',  'V45A'),
    'MUT_49SER':  (54,  'FR2',  'A54S'),
    'MUT_50VAL':  (55,  'FR2',  'I55V'),
    'MUT_74SER':  (83,  'FR3',  'A83S'),
    'MUT_86ARG':  (95,  'FR3',  'K95R'),
    'MUT_87ALA':  (96,  'FR3',  'P96A'),
    'MUT_92VAL':  (101, 'FR3',  'M101V'),
    'MUT_120LEU': (123, 'FR4',  'Q123L'),
    # Set2: I55 saturation
    'MUT_50GLY':  (55,  'FR2',  'I55G'),
    'MUT_50ALA':  (55,  'FR2',  'I55A'),
    'MUT_50SER':  (55,  'FR2',  'I55S'),
    'MUT_50TYR':  (55,  'FR2',  'I55Y'),
    'MUT_50ARG':  (55,  'FR2',  'I55R'),
    'MUT_50PHE':  (55,  'FR2',  'I55F'),
    'MUT_50LEU':  (55,  'FR2',  'I55L'),
    'MUT_50GLN':  (55,  'FR2',  'I55Q'),
    # Set2: P96 saturation
    'MUT_87VAL':  (96,  'FR3',  'P96V'),
    'MUT_87ASP':  (96,  'FR3',  'P96D'),
}


def parse_predictions(filepath):
    """Parse ddG predictions file → list of (mut_name, ddg, std)."""
    results = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ddG:"):
                continue
            parts = line.split()
            name = parts[1]
            if name == "WT" or name == "description":
                continue
            ddg = float(parts[2])
            std = float(parts[3]) if len(parts) > 3 else 0.0
            results.append((name, ddg, std))
    return results


def get_bar_color(ddg):
    """Color by stability classification."""
    if ddg > 2.0:
        return '#d32f2f'   # strong destabilizing - red
    elif ddg > 1.0:
        return '#f57c00'   # moderate - orange
    elif ddg > 0.5:
        return '#fbc02d'   # mild - yellow
    elif ddg > -0.5:
        return '#757575'   # neutral - gray
    else:
        return '#388e3c'   # stabilizing - green


def plot_set1(ax, data):
    """Plot Set1: 13 humanization mutations sorted by IMGT position."""
    # Sort by IMGT position
    sorted_data = sorted(data, key=lambda x: MUTATION_INFO.get(x[0], (999,))[0])

    labels = []
    ddgs = []
    stds = []
    colors = []
    region_labels = []

    for name, ddg, std in sorted_data:
        info = MUTATION_INFO.get(name)
        if info is None:
            continue
        imgt, region, label = info
        labels.append(f"{label}\n(IMGT {imgt})")
        ddgs.append(ddg)
        stds.append(std)
        colors.append(get_bar_color(ddg))
        region_labels.append(region)

    x = np.arange(len(labels))
    bars = ax.bar(x, ddgs, yerr=stds, capsize=3, color=colors, edgecolor='black',
                  linewidth=0.5, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, ha='center')
    ax.set_ylabel('ΔΔG (REU)', fontsize=12)
    ax.set_title('(A) Set1: 13 Humanization Mutations', fontsize=13, fontweight='bold')
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.axhline(y=2.0, color='red', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.axhline(y=-0.5, color='green', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.grid(axis='y', alpha=0.3, zorder=0)

    # Region background shading
    prev_region = None
    region_starts = []
    for i, r in enumerate(region_labels):
        if r != prev_region:
            region_starts.append((i, r))
            prev_region = r

    region_colors = {'FR1': '#fff3e0', 'FR2': '#e3f2fd', 'FR3': '#f3e5f5', 'FR4': '#e8f5e9',
                     'CDR1': '#ffebee', 'CDR2': '#ffebee', 'CDR3': '#ffebee'}
    for idx, (start, region) in enumerate(region_starts):
        end = region_starts[idx + 1][0] if idx + 1 < len(region_starts) else len(labels)
        ax.axvspan(start - 0.5, end - 0.5, alpha=0.15, color=region_colors.get(region, '#f5f5f5'), zorder=0)
        mid = (start + end - 1) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.95, region, ha='center', fontsize=8,
                fontweight='bold', color='gray')


def plot_set2(ax, data):
    """Plot Set2: I55 + P96 saturation, split into two groups."""
    i55_data = [(n, d, s) for n, d, s in data if n.startswith('MUT_50')]
    p96_data = [(n, d, s) for n, d, s in data if n.startswith('MUT_87')]

    # Sort by ddG
    i55_data.sort(key=lambda x: x[1])
    p96_data.sort(key=lambda x: x[1])

    all_data = i55_data + p96_data
    labels = []
    ddgs = []
    stds = []
    colors = []

    for name, ddg, std in all_data:
        info = MUTATION_INFO.get(name)
        label = info[2] if info else name
        labels.append(label)
        ddgs.append(ddg)
        stds.append(std)
        colors.append(get_bar_color(ddg))

    x = np.arange(len(labels))
    ax.bar(x, ddgs, yerr=stds, capsize=3, color=colors, edgecolor='black',
           linewidth=0.5, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9, ha='center')
    ax.set_ylabel('ΔΔG (REU)', fontsize=12)
    ax.set_title('(B) Set2: Position Saturation (IMGT 55, 96)', fontsize=13, fontweight='bold')
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.axhline(y=2.0, color='red', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.axhline(y=-0.5, color='green', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.grid(axis='y', alpha=0.3, zorder=0)

    # Divider between I55 and P96
    if i55_data and p96_data:
        ax.axvline(x=len(i55_data) - 0.5, color='gray', linewidth=1, linestyle=':')
        ax.text(len(i55_data) / 2 - 0.5, ax.get_ylim()[1] * 0.95,
                'IMGT 55 (I→X)', ha='center', fontsize=9, fontweight='bold', color='gray')
        ax.text(len(i55_data) + len(p96_data) / 2 - 0.5, ax.get_ylim()[1] * 0.95,
                'IMGT 96 (P→X)', ha='center', fontsize=9, fontweight='bold', color='gray')


def main():
    set1 = parse_predictions(SET1_FILE)
    set2 = parse_predictions(SET2_FILE)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [1, 0.8]})

    plot_set1(ax1, set1)
    plot_set2(ax2, set2)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d32f2f', edgecolor='black', label='Strongly destabilizing (>2.0)'),
        Patch(facecolor='#f57c00', edgecolor='black', label='Moderately destabilizing (1.0–2.0)'),
        Patch(facecolor='#fbc02d', edgecolor='black', label='Mildly destabilizing (0.5–1.0)'),
        Patch(facecolor='#757575', edgecolor='black', label='Neutral (−0.5–0.5)'),
        Patch(facecolor='#388e3c', edgecolor='black', label='Stabilizing (<−0.5)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    out_png = os.path.join(FIGURES_DIR, "04_ddg_barplot.png")
    out_pdf = os.path.join(FIGURES_DIR, "04_ddg_barplot.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")

    # Save plot data as CSV
    out_csv = os.path.join(SCRIPT_DIR, "04_ddg_plot_data.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["set", "rosetta_name", "label", "imgt_pos", "region", "ddg", "std"])
        for name, ddg, std in set1:
            info = MUTATION_INFO.get(name, (None, None, name))
            writer.writerow(["set1", name, info[2], info[0], info[1], f"{ddg:.3f}", f"{std:.3f}"])
        for name, ddg, std in set2:
            info = MUTATION_INFO.get(name, (None, None, name))
            writer.writerow(["set2", name, info[2], info[0], info[1], f"{ddg:.3f}", f"{std:.3f}"])
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
