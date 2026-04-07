#!/usr/bin/env python3
"""
Destabilizing mutation position의 VH3 germline AA 분포 stacked bar plot

Set1에서 destabilizing으로 판정된 6개 position (IMGT 1, 15, 45, 54, 55, 96)에 대해
human VH3 germline에서 해당 위치의 아미노산 비율을 stacked bar plot으로 시각화.
이 분포가 Set2 대안 AA 선택 및 최종 reject/accept 판단의 근거가 됨.

Input:
    - ../03_vh3_family_analysis/04_VH3_imgt_aa_proportion.csv : VH3 AA proportions
Output:
    - figures/06_vh3_proportion_destabilizing.png : stacked bar plot (PNG)
    - figures/06_vh3_proportion_destabilizing.pdf : stacked bar plot (PDF)
    - 06_vh3_proportion_destabilizing.csv         : plot data (CSV)

Usage:
    python 06_vh3_proportion_destabilizing.py
"""

import csv
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.pad': 4,
    'ytick.major.pad': 4,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

VH3_FILE = os.path.join(SCRIPT_DIR, "..", "03_vh3_family_analysis", "04_VH3_imgt_aa_proportion.csv")

# 6 destabilizing positions: (IMGT, nanobody_WT, humanization_target, ddG, region)
POSITIONS = [
    (1,  'Q', 'E', +2.09, 'FR1'),
    (15, 'A', 'P', +6.04, 'FR1'),
    (45, 'V', 'A', +2.76, 'FR2'),
    (54, 'A', 'S', +4.16, 'FR2'),
    (55, 'I', 'V', +2.52, 'FR2'),
    (96, 'P', 'A', +3.44, 'FR3'),
]

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# Muted, colorblind-friendly palette (Tableau 20 inspired)
AA_COLORS = {
    'A': '#4e79a7', 'C': '#59a14f', 'D': '#e15759', 'E': '#76b7b2',
    'F': '#f28e2b', 'G': '#b07aa1', 'H': '#ff9da7', 'I': '#9c755f',
    'K': '#bab0ac', 'L': '#edc948', 'M': '#8cd17d', 'N': '#c49c94',
    'P': '#d37295', 'Q': '#a0cbe8', 'R': '#b6992d', 'S': '#86bcb6',
    'T': '#499894', 'V': '#f1ce63', 'W': '#79706e', 'Y': '#d4a6c8',
}

REGION_COLORS = {
    'FR1': '#3f51b5', 'FR2': '#00897b', 'FR3': '#8e24aa',
}

# Threshold below which labels are placed outside the bar segment
LABEL_MIN_PROPORTION = 0.06
LABEL_SHOW_PROPORTION = 0.03


def load_vh3_proportions(filepath):
    """Load VH3 AA proportions -> {position_str: {AA: proportion}}."""
    data = {}
    with open(filepath) as f:
        for row in csv.DictReader(f):
            pos = row['position']
            props = {}
            for aa in AA_ORDER:
                val = float(row.get(aa, 0))
                if val > 0:
                    props[aa] = val
            data[pos] = props
    return data


def label_color_for_aa(aa):
    """Return a contrasting text color for the given amino acid segment."""
    dark_bg_aas = {'W', 'R', 'I', 'A', 'T', 'D'}
    if aa in dark_bg_aas:
        return 'white'
    return '#1a1a1a'


def main():
    vh3 = load_vh3_proportions(VH3_FILE)

    fig, ax = plt.subplots(figsize=(10, 6.2))

    x_positions = np.arange(len(POSITIONS))
    bar_width = 0.52
    all_plot_data = []

    for idx, (imgt, wt, target, ddg, region) in enumerate(POSITIONS):
        props = vh3.get(str(imgt), {})
        sorted_aas = sorted(props.items(), key=lambda kv: -kv[1])

        bottom = 0.0
        for aa, prop in sorted_aas:
            color = AA_COLORS.get(aa, '#cccccc')

            is_wt = (aa == wt)
            is_target = (aa == target)
            hatch = '////' if is_wt else None
            edgecolor = '#222222' if is_wt or is_target else 'white'
            linewidth = 2.0 if is_wt or is_target else 0.5

            ax.bar(idx, prop, bottom=bottom, color=color, edgecolor=edgecolor,
                   linewidth=linewidth, width=bar_width, hatch=hatch,
                   zorder=2)

            if prop >= LABEL_SHOW_PROPORTION:
                label_text = f"{aa} {prop:.0%}"
                if prop >= LABEL_MIN_PROPORTION:
                    ax.text(idx, bottom + prop / 2, label_text,
                            ha='center', va='center',
                            fontsize=8.5 if prop >= 0.15 else 7,
                            fontweight='bold',
                            color=label_color_for_aa(aa), zorder=3)
                else:
                    # Small segments: place label to the right of the bar
                    ax.annotate(
                        label_text,
                        xy=(idx + bar_width / 2, bottom + prop / 2),
                        xytext=(6, 0), textcoords='offset points',
                        ha='left', va='center', fontsize=6.5,
                        fontweight='bold', color='#444444',
                        arrowprops=dict(arrowstyle='-', color='#aaaaaa',
                                        lw=0.5, shrinkA=0, shrinkB=2),
                        zorder=3)

            bottom += prop

            all_plot_data.append({
                'imgt_pos': imgt, 'region': region, 'nanobody_wt': wt,
                'humanization_target': target, 'set1_ddg': ddg,
                'amino_acid': aa, 'vh3_proportion': prop,
            })

    # X-axis labels
    ax.set_xticks(x_positions)
    x_labels = [
        f"IMGT {imgt}\n{wt}$\\rightarrow${target}  (+{ddg:.1f})"
        for imgt, wt, target, ddg, region in POSITIONS
    ]
    ax.set_xticklabels(x_labels, fontsize=9, fontweight='bold')

    # Region badges below x-axis
    for idx, (_, _, _, _, region) in enumerate(POSITIONS):
        ax.annotate(region, xy=(idx, -0.07), xycoords=('data', 'axes fraction'),
                    ha='center', va='top', fontsize=7.5, fontweight='bold',
                    color='white', bbox=dict(
                        boxstyle='round,pad=0.3', facecolor=REGION_COLORS[region],
                        edgecolor='none', alpha=0.9))

    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.5, len(POSITIONS) - 0.5)
    ax.set_ylabel('VH3 Germline AA Proportion', fontsize=11, fontweight='bold',
                  labelpad=8)
    ax.set_title('VH3 Germline AA Distribution at Destabilizing Positions',
                 fontsize=13, fontweight='bold', pad=14)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#666666')
    ax.spines['bottom'].set_color('#666666')
    ax.tick_params(colors='#444444')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0%}'))

    ax.yaxis.grid(True, linestyle='--', alpha=0.25, color='#888888', zorder=0)
    ax.set_axisbelow(True)

    legend_elements = [
        mpatches.Patch(facecolor='#dddddd', edgecolor='#222222', hatch='////',
                       linewidth=1.5, label='Nanobody WT residue'),
        mpatches.Patch(facecolor='#dddddd', edgecolor='#222222',
                       linewidth=1.5, label='Humanization target'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8.5,
              framealpha=0.92, edgecolor='#cccccc', fancybox=True,
              borderpad=0.8, handlelength=1.8)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)

    out_png = os.path.join(FIGURES_DIR, "06_vh3_proportion_destabilizing.png")
    out_pdf = os.path.join(FIGURES_DIR, "06_vh3_proportion_destabilizing.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    plt.close(fig)

    out_csv = os.path.join(SCRIPT_DIR, "06_vh3_proportion_destabilizing.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            'imgt_pos', 'region', 'nanobody_wt', 'humanization_target',
            'set1_ddg', 'amino_acid', 'vh3_proportion',
        ])
        writer.writeheader()
        writer.writerows(all_plot_data)
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
