#!/usr/bin/env python3
"""
인간화 전후 nanobody의 humanness score 계산 및 시각화

두 가지 지표로 인간화 효과를 정량화:
  1. FR Identity to IGHV3-66 : FR 잔기 중 IGHV3-66*01 germline과 일치하는 비율
  2. Position-Weighted Humanness (PWH) : 각 FR position에서 해당 AA의 VH3 family 비율 평균

원본(WT) nanobody와 인간화(Humanized) nanobody를 비교하여
position별 humanness 변화를 시각화한다.

Input:
    - ../02_rosetta_ddg/01_imgt_pdb_mapping.csv           : IMGT ↔ PDB mapping (WT)
    - ../03_vh3_family_analysis/04_VH3_imgt_aa_proportion.csv : VH3 AA proportions
    - ../03_vh3_family_analysis/03_VH3_unique_imgt.csv    : IGHV3-66 IMGT residues
Output:
    - 01_humanness_score.csv              : position별 humanness 데이터
    - 01_humanness_summary.txt            : 요약 보고서
    - figures/01_humanness_comparison.png  : humanness score 비교 bar plot
    - figures/01_humanness_comparison.pdf  : PDF

Usage:
    python 01_humanness_score.py
"""

import csv
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

MAPPING_FILE = os.path.join(SCRIPT_DIR, "..", "02_rosetta_ddg", "01_imgt_pdb_mapping.csv")
VH3_PROP_FILE = os.path.join(SCRIPT_DIR, "..", "03_vh3_family_analysis",
                             "04_VH3_imgt_aa_proportion.csv")
VH3_IMGT_FILE = os.path.join(SCRIPT_DIR, "..", "03_vh3_family_analysis",
                             "03_VH3_unique_imgt.csv")

# Original nanobody sequence
WT_SEQ = "QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS"

# 8 accepted mutations: (IMGT, PDB_pos, WT_aa, MUT_aa)
ACCEPTED_MUTATIONS = [
    (5,   5,   'Q', 'V'),
    (12,  11,  'S', 'L'),
    (40,  35,  'G', 'S'),
    (55,  50,  'I', 'L'),
    (83,  74,  'A', 'S'),
    (95,  86,  'K', 'R'),
    (101, 92,  'M', 'V'),
    (123, 120, 'Q', 'L'),
]

FR_REGIONS = {'FR1', 'FR2', 'FR3', 'FR4'}

REGION_COLORS = {
    'FR1': '#5c6bc0', 'FR2': '#26a69a', 'FR3': '#ab47bc', 'FR4': '#66bb6a',
}


def build_humanized_seq():
    """Apply 8 accepted mutations to WT sequence."""
    seq = list(WT_SEQ)
    for imgt, pdb_pos, wt, mut in ACCEPTED_MUTATIONS:
        assert seq[pdb_pos - 1] == wt, f"Mismatch at PDB {pdb_pos}: expected {wt}, got {seq[pdb_pos - 1]}"
        seq[pdb_pos - 1] = mut
    return ''.join(seq)


def load_mapping():
    """Load IMGT mapping for FR positions."""
    positions = []
    with open(MAPPING_FILE) as f:
        for row in csv.DictReader(f):
            positions.append({
                'pdb_pos': int(row['pdb_pos']),
                'aa': row['aa'],
                'imgt_pos': row['imgt_pos'],
                'insertion': row['insertion'],
                'region': row['region'],
            })
    return positions


def load_vh3_proportions():
    """Load VH3 AA proportions → {imgt_pos_str: {AA: proportion}}."""
    data = {}
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    with open(VH3_PROP_FILE) as f:
        for row in csv.DictReader(f):
            pos = row['position']
            props = {}
            for aa in aa_list:
                val = float(row.get(aa, 0))
                if val > 0:
                    props[aa] = val
            data[pos] = props
    return data


def load_ighv366_residues():
    """Load IGHV3-66*01 IMGT-numbered residues."""
    residues = {}
    with open(VH3_IMGT_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'IGHV3-66*01' in row.get('all_alleles', ''):
                for key, val in row.items():
                    if key[0].isdigit() and val not in ('-', ''):
                        residues[key] = val
                break
    return residues


def main():
    humanized_seq = build_humanized_seq()
    mapping = load_mapping()
    vh3_props = load_vh3_proportions()
    ighv366 = load_ighv366_residues()

    print(f"WT seq:        {WT_SEQ}")
    print(f"Humanized seq: {humanized_seq}")
    print(f"Differences:   {sum(1 for a, b in zip(WT_SEQ, humanized_seq) if a != b)} positions")
    print()

    # Build per-position data
    results = []
    for m in mapping:
        region = m['region']
        if region not in FR_REGIONS:
            continue
        if m['insertion']:
            continue

        pdb_pos = m['pdb_pos']
        imgt = m['imgt_pos']
        wt_aa = WT_SEQ[pdb_pos - 1]
        hum_aa = humanized_seq[pdb_pos - 1]

        # IGHV3-66 residue at this position
        germline_aa = ighv366.get(imgt, '')

        # VH3 family proportion
        pos_props = vh3_props.get(imgt, {})
        wt_vh3_prop = pos_props.get(wt_aa, 0.0)
        hum_vh3_prop = pos_props.get(hum_aa, 0.0)

        # Identity to IGHV3-66
        wt_match = 1 if (germline_aa and wt_aa == germline_aa) else 0
        hum_match = 1 if (germline_aa and hum_aa == germline_aa) else 0

        mutated = wt_aa != hum_aa

        results.append({
            'pdb_pos': pdb_pos,
            'imgt_pos': imgt,
            'region': region,
            'wt_aa': wt_aa,
            'hum_aa': hum_aa,
            'germline_aa': germline_aa,
            'wt_germline_match': wt_match,
            'hum_germline_match': hum_match,
            'wt_vh3_prop': wt_vh3_prop,
            'hum_vh3_prop': hum_vh3_prop,
            'mutated': mutated,
        })

    # --- Metrics ---
    n_fr = len(results)
    n_with_germline = sum(1 for r in results if r['germline_aa'])

    # 1. FR Identity to IGHV3-66
    wt_identity = sum(r['wt_germline_match'] for r in results) / n_with_germline
    hum_identity = sum(r['hum_germline_match'] for r in results) / n_with_germline

    # 2. Position-Weighted Humanness (mean VH3 proportion across FR)
    wt_pwh = np.mean([r['wt_vh3_prop'] for r in results])
    hum_pwh = np.mean([r['hum_vh3_prop'] for r in results])

    print(f"=== Humanness Score ===")
    print(f"FR positions evaluated: {n_fr} (with germline data: {n_with_germline})")
    print()
    print(f"1. FR Identity to IGHV3-66*01:")
    print(f"   WT:        {wt_identity:.1%} ({sum(r['wt_germline_match'] for r in results)}/{n_with_germline})")
    print(f"   Humanized: {hum_identity:.1%} ({sum(r['hum_germline_match'] for r in results)}/{n_with_germline})")
    print(f"   Change:    +{(hum_identity - wt_identity):.1%}")
    print()
    print(f"2. Position-Weighted Humanness (VH3 family):")
    print(f"   WT:        {wt_pwh:.4f}")
    print(f"   Humanized: {hum_pwh:.4f}")
    print(f"   Change:    {hum_pwh - wt_pwh:+.4f}")

    # --- Save CSV ---
    out_csv = os.path.join(SCRIPT_DIR, "01_humanness_score.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pdb_pos', 'imgt_pos', 'region', 'wt_aa', 'hum_aa', 'germline_aa',
            'wt_germline_match', 'hum_germline_match', 'wt_vh3_prop', 'hum_vh3_prop', 'mutated',
        ])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {out_csv}")

    # --- Save summary text ---
    out_txt = os.path.join(SCRIPT_DIR, "01_humanness_summary.txt")
    with open(out_txt, "w") as f:
        f.write("Nanobody Humanness Score Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"WT sequence:        {WT_SEQ}\n")
        f.write(f"Humanized sequence: {humanized_seq}\n")
        f.write(f"Mutations applied:  {sum(1 for r in results if r['mutated'])} / {n_fr} FR positions\n\n")
        f.write(f"1. FR Identity to IGHV3-66*01\n")
        f.write(f"   WT:        {wt_identity:.1%} ({sum(r['wt_germline_match'] for r in results)}/{n_with_germline})\n")
        f.write(f"   Humanized: {hum_identity:.1%} ({sum(r['hum_germline_match'] for r in results)}/{n_with_germline})\n")
        f.write(f"   Change:    +{(hum_identity - wt_identity):.1%}\n\n")
        f.write(f"2. Position-Weighted Humanness (VH3 family mean proportion)\n")
        f.write(f"   WT:        {wt_pwh:.4f}\n")
        f.write(f"   Humanized: {hum_pwh:.4f}\n")
        f.write(f"   Change:    {hum_pwh - wt_pwh:+.4f}\n\n")
        f.write(f"Accepted mutations:\n")
        for r in results:
            if r['mutated']:
                f.write(f"  IMGT {r['imgt_pos']:>3} ({r['region']}): "
                        f"{r['wt_aa']}→{r['hum_aa']}  "
                        f"VH3 prop: {r['wt_vh3_prop']:.3f}→{r['hum_vh3_prop']:.3f}  "
                        f"germline({r['germline_aa']}): {'match' if r['hum_germline_match'] else 'no match'}\n")
    print(f"Saved: {out_txt}")

    # --- Visualization: vertical grouped bar chart ---
    fig, ax = plt.subplots(figsize=(6, 5))

    metrics = ['FR Identity\n(IGHV3-66)', 'Position-Weighted\nHumanness (VH3)']
    wt_vals = [wt_identity * 100, wt_pwh * 100]
    hum_vals = [hum_identity * 100, hum_pwh * 100]

    x = np.arange(len(metrics))
    width = 0.32

    bars_wt = ax.bar(x - width/2, wt_vals, width, label='WT (original)',
                      color='#ef5350', alpha=0.85, edgecolor='white', linewidth=0.5, zorder=2)
    bars_hum = ax.bar(x + width/2, hum_vals, width, label='Humanized',
                       color='#42a5f5', alpha=0.85, edgecolor='white', linewidth=0.5, zorder=2)

    for bar, val in zip(bars_wt, wt_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold',
                color='#ef5350')
    for bar, val in zip(bars_hum, hum_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold',
                color='#42a5f5')

    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=11, fontweight='bold')
    ax.set_ylabel('Score (%)', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.set_title('Humanness Score: WT vs Humanized Nanobody',
                 fontsize=13, fontweight='bold', pad=12)
    ax.legend(fontsize=10, loc='upper left', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.2, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()

    out_png = os.path.join(FIGURES_DIR, "01_humanness_comparison.png")
    out_pdf = os.path.join(FIGURES_DIR, "01_humanness_comparison.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    plt.close(fig)


if __name__ == "__main__":
    main()
