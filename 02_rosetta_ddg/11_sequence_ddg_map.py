#!/usr/bin/env python3
"""
나노바디 서열 위에 ddG 값을 매핑한 lollipop plot 생성

Visualizes ddG values mapped onto the nanobody sequence with FR/CDR
region annotation. Each humanization candidate position is shown as a
lollipop (stem + circle), colored by stability classification. The
sequence is displayed below with region shading.

Input:
    - 02_ddg_predictions_final.out : Set1 ddG predictions
    - 01_imgt_pdb_mapping.csv      : IMGT ↔ PDB position mapping
Output:
    - figures/11_sequence_ddg_map.png : lollipop plot (PNG)
    - figures/11_sequence_ddg_map.pdf : lollipop plot (PDF)
    - 11_sequence_ddg_map_data.csv    : plot data (CSV)

Usage:
    python 11_sequence_ddg_map.py
"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

DDG_FILE = os.path.join(SCRIPT_DIR, "02_ddg_predictions_final.out")
MAPPING_FILE = os.path.join(SCRIPT_DIR, "01_imgt_pdb_mapping.csv")

TARGET_SEQ = "QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS"

# Rosetta MUT name → (IMGT, label, final_mutation)
# final_mutation: what was actually decided (I55→L instead of I55→V)
MUTATION_INFO = {
    'MUT_1GLU':   (1,   'Q1E',   False),
    'MUT_5VAL':   (5,   'Q5V',   True),
    'MUT_11LEU':  (12,  'S12L',  True),
    'MUT_14PRO':  (15,  'A15P',  False),
    'MUT_35SER':  (40,  'G40S',  True),
    'MUT_40ALA':  (45,  'V45A',  False),
    'MUT_49SER':  (54,  'A54S',  False),
    'MUT_50VAL':  (55,  'I55V→L', True),  # replaced by I55L
    'MUT_74SER':  (83,  'A83S',  True),
    'MUT_86ARG':  (95,  'K95R',  True),
    'MUT_87ALA':  (96,  'P96A',  False),
    'MUT_92VAL':  (101, 'M101V', True),
    'MUT_120LEU': (123, 'Q123L', True),
}

REGION_COLORS = {
    'FR1': '#e8eaf6', 'CDR1': '#ffcdd2', 'FR2': '#e0f2f1',
    'CDR2': '#fff9c4', 'FR3': '#f3e5f5', 'CDR3': '#ffccbc', 'FR4': '#e8f5e9',
}


def parse_predictions(filepath):
    results = {}
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ddG:"):
                continue
            parts = line.split()
            name = parts[1]
            if name in ("WT", "description"):
                continue
            results[name] = (float(parts[2]), float(parts[3]) if len(parts) > 3 else 0.0)
    return results


def load_mapping(filepath):
    """Load IMGT-PDB mapping → list of (pdb_pos, aa, imgt_pos, region)."""
    mapping = []
    with open(filepath) as f:
        for row in csv.DictReader(f):
            mapping.append({
                'pdb_pos': int(row['pdb_pos']),
                'aa': row['aa'],
                'imgt_pos': row['imgt_pos'],
                'insertion': row['insertion'],
                'region': row['region'],
            })
    return mapping


def get_bar_color(ddg):
    if ddg > 2.0:
        return '#d32f2f'
    elif ddg > 1.0:
        return '#f57c00'
    elif ddg > 0.5:
        return '#fbc02d'
    elif ddg > -0.5:
        return '#757575'
    else:
        return '#388e3c'


def main():
    ddg_data = parse_predictions(DDG_FILE)
    mapping = load_mapping(MAPPING_FILE)

    # Build PDB pos → (IMGT, region) lookup
    pdb_to_info = {}
    for m in mapping:
        pdb_to_info[m['pdb_pos']] = m

    # Build mutation positions: pdb_pos → (ddg, std, label, accepted)
    mut_positions = {}
    for rosetta_name, (ddg, std) in ddg_data.items():
        info = MUTATION_INFO.get(rosetta_name)
        if info is None:
            continue
        imgt, label, accepted = info
        # Find PDB pos for this IMGT
        for m in mapping:
            if not m['insertion'] and m['imgt_pos'] == str(imgt):
                mut_positions[m['pdb_pos']] = (ddg, std, label, accepted, imgt, m['region'])
                break

    # --- Figure ---
    fig, (ax_lollipop, ax_seq) = plt.subplots(
        2, 1, figsize=(20, 6), gridspec_kw={'height_ratios': [3, 1]},
        sharex=True
    )

    seq_len = len(TARGET_SEQ)
    x_all = np.arange(1, seq_len + 1)

    # Region shading on both axes
    prev_region = None
    region_start = 1
    for m in mapping:
        region = m['region']
        pos = m['pdb_pos']
        if region != prev_region:
            if prev_region is not None:
                for ax in (ax_lollipop, ax_seq):
                    ax.axvspan(region_start - 0.5, pos - 0.5, alpha=0.3,
                               color=REGION_COLORS.get(prev_region, '#f5f5f5'), zorder=0)
                # Region label on lollipop
                mid = (region_start + pos - 1) / 2
                ax_lollipop.text(mid, -5.5, prev_region, ha='center', fontsize=7,
                                 fontweight='bold', color='#555555')
            region_start = pos
            prev_region = region
    # Last region
    for ax in (ax_lollipop, ax_seq):
        ax.axvspan(region_start - 0.5, seq_len + 0.5, alpha=0.3,
                   color=REGION_COLORS.get(prev_region, '#f5f5f5'), zorder=0)
    mid = (region_start + seq_len) / 2
    ax_lollipop.text(mid, -5.5, prev_region, ha='center', fontsize=7,
                     fontweight='bold', color='#555555')

    # Lollipop plot
    ax_lollipop.axhline(y=0, color='black', linewidth=0.8, zorder=1)
    ax_lollipop.axhline(y=2.0, color='red', linewidth=0.6, linestyle='--', alpha=0.4)
    ax_lollipop.axhline(y=-0.5, color='green', linewidth=0.6, linestyle='--', alpha=0.4)

    for pdb_pos, (ddg, std, label, accepted, imgt, region) in mut_positions.items():
        color = get_bar_color(ddg)
        marker = 'o' if accepted else 'X'
        markersize = 10 if accepted else 9

        # Stem
        ax_lollipop.plot([pdb_pos, pdb_pos], [0, ddg], color=color, linewidth=2, zorder=2)
        # Circle/X
        ax_lollipop.plot(pdb_pos, ddg, marker=marker, color=color, markersize=markersize,
                         markeredgecolor='black', markeredgewidth=0.8, zorder=3)
        # Error bar
        ax_lollipop.plot([pdb_pos, pdb_pos], [ddg - std, ddg + std],
                         color=color, linewidth=1, alpha=0.5, zorder=2)
        # Label
        va = 'bottom' if ddg >= 0 else 'top'
        offset = 0.4 if ddg >= 0 else -0.4
        ax_lollipop.text(pdb_pos, ddg + offset + (std if ddg >= 0 else -std),
                         label, ha='center', va=va, fontsize=7, fontweight='bold',
                         rotation=45)

    ax_lollipop.set_ylabel('ΔΔG (REU)', fontsize=12)
    ax_lollipop.set_title('Humanization Mutation Stability Map on Nanobody Sequence',
                          fontsize=14, fontweight='bold')
    ax_lollipop.set_xlim(0.5, seq_len + 0.5)

    # Sequence display
    for i, aa in enumerate(TARGET_SEQ):
        pos = i + 1
        fontweight = 'bold' if pos in mut_positions else 'normal'
        color = 'red' if (pos in mut_positions and not mut_positions[pos][3]) else \
                'blue' if pos in mut_positions else 'black'
        ax_seq.text(pos, 0.5, aa, ha='center', va='center', fontsize=5.5,
                    fontweight=fontweight, color=color, fontfamily='monospace')

    # Position numbers every 10
    for i in range(0, seq_len, 10):
        ax_seq.text(i + 1, -0.1, str(i + 1), ha='center', va='top', fontsize=6, color='gray')

    ax_seq.set_ylim(-0.5, 1.0)
    ax_seq.set_yticks([])
    ax_seq.set_xlabel('PDB Sequential Position', fontsize=11)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#d32f2f', label='Destabilizing (>2.0)'),
        mpatches.Patch(facecolor='#388e3c', label='Stabilizing (<−0.5)'),
        mpatches.Patch(facecolor='#757575', label='Neutral (−0.5–0.5)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=8, markeredgecolor='black', label='Accepted'),
        plt.Line2D([0], [0], marker='X', color='w', markerfacecolor='gray',
                   markersize=8, markeredgecolor='black', label='Rejected'),
    ]
    ax_lollipop.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.9)

    plt.tight_layout()

    out_png = os.path.join(FIGURES_DIR, "11_sequence_ddg_map.png")
    out_pdf = os.path.join(FIGURES_DIR, "11_sequence_ddg_map.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")

    # Save plot data CSV
    out_csv = os.path.join(SCRIPT_DIR, "11_sequence_ddg_map_data.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "aa", "imgt_pos", "region", "mutation", "ddg", "std", "accepted"])
        for pdb_pos in sorted(mut_positions.keys()):
            ddg, std, label, accepted, imgt, region = mut_positions[pdb_pos]
            writer.writerow([pdb_pos, TARGET_SEQ[pdb_pos - 1], imgt, region,
                             label, f"{ddg:.3f}", f"{std:.3f}", accepted])
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
