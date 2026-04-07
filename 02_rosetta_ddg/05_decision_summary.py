#!/usr/bin/env python3
"""
ddG 결과와 VH3 germline AA proportion을 통합하여 최종 인간화 결정 요약 테이블 생성

Combines Rosetta Cartesian ddG predictions with human VH3 family amino acid
proportions at each IMGT position to produce a unified decision table showing
the rationale for accepting or rejecting each humanization mutation.

Input:
    - 02_ddg_predictions_final.out : Set1 ddG predictions (13 mutations)
    - 02_ddg_predictions_set2.out  : Set2 ddG predictions (I55/P96 saturation)
    - ../03_vh3_family_analysis/VH3_imgt_aa_proportion.csv : VH3 AA proportions
Output:
    - 05_humanization_decision.csv : final decision table
    - 05_humanization_decision.txt : formatted text report
    - figures/05_decision_table.png : visual summary (PNG)
    - figures/05_decision_table.pdf : visual summary (PDF)

Usage:
    python 05_decision_summary.py
"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

SET1_FILE = os.path.join(SCRIPT_DIR, "02_ddg_predictions_final.out")
SET2_FILE = os.path.join(SCRIPT_DIR, "02_ddg_predictions_set2.out")
VH3_FILE = os.path.join(SCRIPT_DIR, "..", "03_vh3_family_analysis", "VH3_imgt_aa_proportion.csv")

# 13 humanization mutations: (IMGT, WT_1letter, MUT_1letter, rosetta_key, region)
MUTATIONS = [
    (1,   'Q', 'E', 'MUT_1GLU',   'FR1'),
    (5,   'Q', 'V', 'MUT_5VAL',   'FR1'),
    (12,  'S', 'L', 'MUT_11LEU',  'FR1'),
    (15,  'A', 'P', 'MUT_14PRO',  'FR1'),
    (40,  'G', 'S', 'MUT_35SER',  'FR2'),
    (45,  'V', 'A', 'MUT_40ALA',  'FR2'),
    (54,  'A', 'S', 'MUT_49SER',  'FR2'),
    (55,  'I', 'V', 'MUT_50VAL',  'FR2'),
    (83,  'A', 'S', 'MUT_74SER',  'FR3'),
    (95,  'K', 'R', 'MUT_86ARG',  'FR3'),
    (96,  'P', 'A', 'MUT_87ALA',  'FR3'),
    (101, 'M', 'V', 'MUT_92VAL',  'FR3'),
    (123, 'Q', 'L', 'MUT_120LEU', 'FR4'),
]

# Final decisions (hardcoded based on analysis)
# True = accept mutation (humanize), False = reject (keep original)
FINAL_DECISIONS = {
    1:   (False, "Destabilizing; Q present in 14.9% of VH3"),
    5:   (True,  "Neutral ddG"),
    12:  (True,  "Stabilizing ddG"),
    15:  (False, "Strongly destabilizing (backbone rigidity)"),
    40:  (True,  "Stabilizing ddG"),
    45:  (False, "Destabilizing; ddG too high"),
    54:  (False, "Destabilizing; A present in 11.9% of VH3"),
    55:  (True,  "I55L selected from Set2 saturation (ddG = -0.25)"),
    83:  (True,  "Stabilizing ddG"),
    95:  (True,  "Neutral ddG (conservative K→R)"),
    96:  (False, "All variants destabilizing (Set2: P96V +3.68, P96D +5.62)"),
    101: (True,  "Stabilizing ddG"),
    123: (True,  "Stabilizing ddG"),
}


def parse_predictions(filepath):
    """Parse ddG predictions → dict of {rosetta_name: (ddg, std)}."""
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


def load_vh3_proportions(filepath):
    """Load VH3 AA proportions → dict of {imgt_pos: {AA: proportion}}."""
    props = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos = row['position']
            if not pos[0].isdigit():
                continue
            # Only main positions (no insertion letters)
            if any(c.isalpha() for c in pos):
                continue
            pos_int = int(pos)
            aa_props = {}
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                val = float(row.get(aa, 0))
                if val > 0:
                    aa_props[aa] = val
            props[pos_int] = aa_props
    return props


def main():
    set1_ddg = parse_predictions(SET1_FILE)
    set2_ddg = parse_predictions(SET2_FILE)
    vh3_props = load_vh3_proportions(VH3_FILE)

    rows = []
    for imgt, wt, mut, rosetta_key, region in MUTATIONS:
        ddg, std = set1_ddg.get(rosetta_key, (None, None))
        accepted, reason = FINAL_DECISIONS[imgt]

        # VH3 proportion of WT residue
        wt_prop = vh3_props.get(imgt, {}).get(wt, 0.0)
        mut_prop = vh3_props.get(imgt, {}).get(mut, 0.0)

        # For IMGT 55: final mutation is I→L not I→V
        final_mut = mut
        final_ddg = ddg
        if imgt == 55:
            final_mut = 'L'
            final_ddg = set2_ddg.get('MUT_50LEU', (None, None))[0]

        rows.append({
            'imgt': imgt,
            'region': region,
            'wt': wt,
            'candidate_mut': mut,
            'final_mut': final_mut,
            'ddg': ddg,
            'std': std,
            'final_ddg': final_ddg,
            'wt_vh3_pct': wt_prop * 100,
            'mut_vh3_pct': mut_prop * 100,
            'accepted': accepted,
            'reason': reason,
        })

    # Save CSV
    out_csv = os.path.join(SCRIPT_DIR, "05_humanization_decision.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            'imgt', 'region', 'wt', 'candidate_mut', 'final_mut',
            'ddg', 'std', 'final_ddg', 'wt_vh3_pct', 'mut_vh3_pct',
            'accepted', 'reason'])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {out_csv}")

    # Save text report
    out_txt = os.path.join(SCRIPT_DIR, "05_humanization_decision.txt")
    with open(out_txt, "w") as f:
        f.write("=" * 100 + "\n")
        f.write("NANOBODY HUMANIZATION - FINAL DECISION SUMMARY\n")
        f.write("=" * 100 + "\n\n")

        f.write(f"{'IMGT':>5} {'Region':>6} {'WT→MUT':>7} {'ddG(REU)':>10} {'WT in VH3':>10} {'MUT in VH3':>11} {'Decision':>10}  Reason\n")
        f.write("-" * 100 + "\n")

        accepted_count = 0
        for r in rows:
            decision = "ACCEPT" if r['accepted'] else "REJECT"
            if r['accepted']:
                accepted_count += 1
            mut_label = r['final_mut'] if r['accepted'] else r['candidate_mut']
            ddg_str = f"{r['final_ddg']:+.2f}" if r['final_ddg'] is not None else "N/A"
            f.write(f"{r['imgt']:>5} {r['region']:>6} {r['wt']}→{mut_label:>1} {ddg_str:>10} "
                    f"{r['wt_vh3_pct']:>9.1f}% {r['mut_vh3_pct']:>10.1f}% "
                    f"{'✓ '+decision if r['accepted'] else '✗ '+decision:>10}  {r['reason']}\n")

        f.write("-" * 100 + "\n")
        f.write(f"\nAccepted: {accepted_count}/13 mutations\n")
        f.write(f"Final humanized mutations: ")
        accepted_muts = [f"{r['wt']}{r['imgt']}{r['final_mut']}" for r in rows if r['accepted']]
        f.write(", ".join(accepted_muts) + "\n")

    print(f"Saved: {out_txt}")

    # Print to console
    with open(out_txt) as f:
        print(f.read())

    # Generate figure: decision table visualization
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.set_axis_off()

    # Table data
    cell_text = []
    cell_colors = []
    for r in rows:
        mut_label = r['final_mut'] if r['accepted'] else r['candidate_mut']
        ddg_str = f"{r['final_ddg']:+.2f}" if r['final_ddg'] is not None else "N/A"
        decision = "✓ Accept" if r['accepted'] else "✗ Reject"

        cell_text.append([
            str(r['imgt']), r['region'], f"{r['wt']}→{mut_label}",
            ddg_str, f"{r['wt_vh3_pct']:.1f}%", f"{r['mut_vh3_pct']:.1f}%",
            decision
        ])

        if r['accepted']:
            cell_colors.append(['#e8f5e9'] * 7)
        else:
            cell_colors.append(['#ffebee'] * 7)

    col_labels = ['IMGT', 'Region', 'Mutation', 'ΔΔG (REU)', 'WT in VH3', 'MUT in VH3', 'Decision']

    table = ax.table(cellText=cell_text, colLabels=col_labels, cellColours=cell_colors,
                     colColours=['#e0e0e0'] * 7, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)

    ax.set_title('Humanization Decision Summary: ddG + VH3 Germline Analysis',
                 fontsize=14, fontweight='bold', pad=20)

    out_png = os.path.join(FIGURES_DIR, "05_decision_table.png")
    out_pdf = os.path.join(FIGURES_DIR, "05_decision_table.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
