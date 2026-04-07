#!/usr/bin/env python3
"""
FR saturation mutagenesis 결과를 분석하고 mutational tolerance landscape 시각화

Parses all position-level ddG results from saturation mutagenesis,
calculates per-position tolerance scores, and generates a heatmap
showing which FR positions are structurally critical vs. tolerant.

Input:
    - saturation/ddg_results/*.ddg : position별 Rosetta ddG raw output
    - 01_imgt_pdb_mapping.csv       : IMGT ↔ PDB mapping
Output:
    - 12_saturation_results.csv         : all ddG values (position × amino acid)
    - 12_position_tolerance.csv         : per-position tolerance summary
    - figures/12_tolerance_heatmap.png  : heatmap (IMGT position × amino acid)
    - figures/12_tolerance_heatmap.pdf  : heatmap (PDF)
    - figures/12_tolerance_barplot.png  : per-position mean |ddG| bar plot
    - figures/12_tolerance_barplot.pdf  : bar plot (PDF)

Usage:
    python 12_analyze_saturation.py
"""

import csv
import glob
import os
import re
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
DDG_DIR = os.path.join(SCRIPT_DIR, "saturation", "ddg_results")
MAPPING_FILE = os.path.join(SCRIPT_DIR, "01_imgt_pdb_mapping.csv")

os.makedirs(FIGURES_DIR, exist_ok=True)

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

REGION_COLORS = {
    'FR1': '#5c6bc0', 'FR2': '#26a69a', 'FR3': '#ab47bc', 'FR4': '#66bb6a',
}


def load_mapping():
    """Load PDB pos → (IMGT, region, WT aa) for FR positions."""
    mapping = {}
    with open(MAPPING_FILE) as f:
        for row in csv.DictReader(f):
            if row['region'] not in ('FR1', 'FR2', 'FR3', 'FR4'):
                continue
            if row['insertion']:
                continue
            mapping[int(row['pdb_pos'])] = {
                'imgt': int(row['imgt_pos']),
                'region': row['region'],
                'wt': row['aa'],
            }
    return mapping


def parse_ddg_file(filepath):
    """Parse a .ddg file → dict of {mut_name: (mean_ddg, std_ddg)}."""
    energies = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            if not line.startswith('COMPLEX:'):
                continue
            parts = line.split()
            mut_name = parts[2].rstrip(':')
            energy = float(parts[3])
            energies[mut_name].append(energy)

    wt_energies = np.array(energies.get('WT', [0]))
    wt_mean = np.mean(wt_energies)

    results = {}
    for name, vals in energies.items():
        if name == 'WT':
            continue
        arr = np.array(vals)
        ddg = np.mean(arr) - wt_mean
        std = np.sqrt(np.std(arr, ddof=1)**2 + np.std(wt_energies, ddof=1)**2) if len(arr) > 1 else 0
        # Extract mutation info from name: MUT_<pos><3letter>
        results[name] = (ddg, std)

    return results


def main():
    mapping = load_mapping()
    ddg_files = sorted(glob.glob(os.path.join(DDG_DIR, "*.ddg")))

    if not ddg_files:
        print(f"No .ddg files found in {DDG_DIR}/")
        print("Run 12_run_saturation_parallel.sh first.")
        print("\nGenerating template output files for demonstration...")
        generate_demo_output(mapping)
        return

    print(f"Found {len(ddg_files)} result files")

    # Parse all results
    # Build matrix: position → {mut_aa: ddg}
    THREE_TO_ONE = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    }

    all_results = []  # (pdb_pos, imgt, region, wt, mut_aa, ddg, std)
    position_data = defaultdict(dict)  # pdb_pos → {mut_aa: ddg}

    for ddg_file in ddg_files:
        results = parse_ddg_file(ddg_file)
        # Extract pdb_pos from filename: mut_FR1_IMGT1_Q1.txt → pdb_pos 1
        basename = os.path.basename(ddg_file).replace('.ddg', '')
        # Parse: mut_<region>_IMGT<imgt>_<wt><pdb_pos>
        match = re.match(r'mut_(\w+)_IMGT(\d+)_(\w)(\d+)', basename)
        if not match:
            continue
        region, imgt_str, wt, pdb_pos_str = match.groups()
        pdb_pos = int(pdb_pos_str)
        imgt = int(imgt_str)

        for mut_name, (ddg, std) in results.items():
            # MUT_<pos><3letter> → extract aa
            mut_match = re.match(r'MUT_\d+(\w{3})', mut_name)
            if not mut_match:
                continue
            mut_3letter = mut_match.group(1).upper()
            # Handle case variations
            for k, v in THREE_TO_ONE.items():
                if k == mut_3letter:
                    mut_aa = v
                    break
            else:
                continue

            all_results.append((pdb_pos, imgt, region, wt, mut_aa, ddg, std))
            position_data[pdb_pos][mut_aa] = ddg

    if not all_results:
        print("No valid results parsed.")
        return

    print(f"Parsed {len(all_results)} mutations across {len(position_data)} positions")

    # Save all results CSV
    out_csv = os.path.join(SCRIPT_DIR, "12_saturation_results.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "imgt_pos", "region", "wt", "mut", "ddg", "std", "mutation_label"])
        for pdb_pos, imgt, region, wt, mut_aa, ddg, std in sorted(all_results):
            label = f"{wt}{imgt}{mut_aa}"
            writer.writerow([pdb_pos, imgt, region, wt, mut_aa, f"{ddg:.3f}", f"{std:.3f}", label])
    print(f"Saved: {out_csv}")

    # Per-position tolerance summary
    tolerance_csv = os.path.join(SCRIPT_DIR, "12_position_tolerance.csv")
    tolerance_data = []
    with open(tolerance_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "imgt_pos", "region", "wt", "mean_abs_ddg",
                         "max_ddg", "min_ddg", "n_tolerant", "n_destabilizing", "tolerance_score"])
        for pdb_pos in sorted(position_data.keys()):
            info = mapping.get(pdb_pos)
            if info is None:
                continue
            ddgs = list(position_data[pdb_pos].values())
            abs_ddgs = [abs(d) for d in ddgs]
            n_tolerant = sum(1 for d in ddgs if abs(d) < 1.0)
            n_destab = sum(1 for d in ddgs if d > 2.0)
            tolerance_score = n_tolerant / len(ddgs) if ddgs else 0

            row = [pdb_pos, info['imgt'], info['region'], info['wt'],
                   f"{np.mean(abs_ddgs):.2f}", f"{max(ddgs):.2f}", f"{min(ddgs):.2f}",
                   n_tolerant, n_destab, f"{tolerance_score:.2f}"]
            writer.writerow(row)
            tolerance_data.append({
                'pdb_pos': pdb_pos, 'imgt': info['imgt'], 'region': info['region'],
                'wt': info['wt'], 'mean_abs_ddg': np.mean(abs_ddgs),
                'tolerance_score': tolerance_score,
            })
    print(f"Saved: {tolerance_csv}")

    # --- Heatmap: IMGT position × amino acid ---
    fr_positions = sorted(position_data.keys())
    imgt_labels = [f"{mapping[p]['wt']}{mapping[p]['imgt']}" for p in fr_positions]

    matrix = np.full((len(AA_ORDER), len(fr_positions)), np.nan)
    for j, pdb_pos in enumerate(fr_positions):
        wt = mapping[pdb_pos]['wt']
        for i, aa in enumerate(AA_ORDER):
            if aa == wt:
                matrix[i, j] = 0.0  # WT = 0
            elif aa in position_data[pdb_pos]:
                matrix[i, j] = position_data[pdb_pos][aa]

    fig, ax = plt.subplots(figsize=(max(20, len(fr_positions) * 0.3), 8))
    cmap = plt.cm.RdYlGn_r
    norm = mcolors.TwoSlopeNorm(vmin=-5, vcenter=0, vmax=10)
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm, interpolation='none')

    ax.set_xticks(range(len(fr_positions)))
    ax.set_xticklabels(imgt_labels, fontsize=5, rotation=90)
    ax.set_yticks(range(len(AA_ORDER)))
    ax.set_yticklabels(AA_ORDER, fontsize=8)
    ax.set_xlabel('FR Position (IMGT)', fontsize=12)
    ax.set_ylabel('Mutant Amino Acid', fontsize=12)
    ax.set_title('Mutational Tolerance Landscape (ΔΔG, REU)', fontsize=14, fontweight='bold')

    plt.colorbar(im, ax=ax, label='ΔΔG (REU)', shrink=0.8)
    plt.tight_layout()

    plt.savefig(os.path.join(FIGURES_DIR, "12_tolerance_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIGURES_DIR, "12_tolerance_heatmap.pdf"), dpi=300, bbox_inches='tight')
    print(f"Saved: figures/12_tolerance_heatmap.png/pdf")
    plt.close()

    # --- Bar plot: per-position mean |ddG| ---
    fig, ax = plt.subplots(figsize=(max(16, len(tolerance_data) * 0.25), 5))

    x = range(len(tolerance_data))
    colors = [REGION_COLORS[t['region']] for t in tolerance_data]
    bars = ax.bar(x, [t['mean_abs_ddg'] for t in tolerance_data], color=colors, edgecolor='black', linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{t['wt']}{t['imgt']}" for t in tolerance_data], fontsize=5, rotation=90)
    ax.set_ylabel('Mean |ΔΔG| (REU)', fontsize=12)
    ax.set_title('Per-Position Mutational Sensitivity (FR regions)', fontsize=14, fontweight='bold')
    ax.axhline(y=2.0, color='red', linestyle='--', linewidth=0.8, alpha=0.5, label='High sensitivity')

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, edgecolor='black', label=r) for r, c in REGION_COLORS.items()]
    ax.legend(handles=legend_elements, fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "12_tolerance_barplot.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIGURES_DIR, "12_tolerance_barplot.pdf"), dpi=300, bbox_inches='tight')
    print(f"Saved: figures/12_tolerance_barplot.png/pdf")
    plt.close()

    # Print top sensitive / tolerant positions
    sorted_tol = sorted(tolerance_data, key=lambda t: t['mean_abs_ddg'], reverse=True)
    print("\nTop 10 most SENSITIVE positions (highest mean |ddG|):")
    for t in sorted_tol[:10]:
        print(f"  IMGT {t['imgt']:>3} ({t['region']}) {t['wt']:>1} — mean |ddG| = {t['mean_abs_ddg']:.2f}")

    print("\nTop 10 most TOLERANT positions (lowest mean |ddG|):")
    for t in sorted_tol[-10:]:
        print(f"  IMGT {t['imgt']:>3} ({t['region']}) {t['wt']:>1} — mean |ddG| = {t['mean_abs_ddg']:.2f}")


def generate_demo_output(mapping):
    """Generate empty template files so the script structure is clear."""
    out_csv = os.path.join(SCRIPT_DIR, "12_saturation_results.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "imgt_pos", "region", "wt", "mut", "ddg", "std", "mutation_label"])
    print(f"Created empty template: {out_csv}")

    tol_csv = os.path.join(SCRIPT_DIR, "12_position_tolerance.csv")
    with open(tol_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "imgt_pos", "region", "wt", "mean_abs_ddg",
                         "max_ddg", "min_ddg", "n_tolerant", "n_destabilizing", "tolerance_score"])
    print(f"Created empty template: {tol_csv}")
    print("\nRun Rosetta saturation first, then re-run this script.")


if __name__ == "__main__":
    main()
