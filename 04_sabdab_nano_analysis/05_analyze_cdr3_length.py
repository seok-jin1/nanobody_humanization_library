#!/usr/bin/env python3
"""
SAbDab 나노바디 구조에서 CDR1/CDR2/CDR3 길이 분포를 계산하고 시각화한다.

Calculate CDR length distributions (CDR1, CDR2, CDR3) from IMGT-numbered
nanobody PDB structures and generate a bar-chart visualization with
descriptive statistics.

Input:
    - sabdab_nano_summary_all.tsv : SAbDab nanobody summary with PDB/chain mapping
    - imgt/*.pdb                  : IMGT-renumbered PDB structure files
Output:
    - 05_cdr_length_distribution.csv        : per-nanobody CDR lengths and sequences
    - figures/05_cdr_length_distribution.png : CDR length distribution bar charts

Usage:
    python 05_analyze_cdr3_length.py
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# --- 설정 ---
SUMMARY_FILE = 'sabdab_nano_summary_all.tsv'
IMGT_DIR = 'imgt'
OUTPUT_CSV = '05_cdr_length_distribution.csv'
OUTPUT_PLOT = 'figures/05_cdr_length_distribution.png'

# IMGT CDR 범위 정의
CDR_REGIONS = {
    'CDR1': (27, 38),
    'CDR2': (56, 65),
    'CDR3': (105, 117),
}

# 3-letter to 1-letter mapping
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X'
}


def load_unique_pairs(summary_file):
    pairs = []
    try:
        df = pd.read_csv(summary_file, sep='\t')
        seen = set()
        for _, row in df.iterrows():
            pdb = str(row['pdb']).lower()
            hchain = str(row['Hchain'])
            if hchain and hchain != 'nan':
                key = (pdb, hchain)
                if key not in seen:
                    seen.add(key)
                    pairs.append(key)
    except Exception as e:
        print(f"Error reading summary file: {e}")
    return pairs


def parse_pdb_residues(pdb_path, target_chain):
    residues = {}
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    if chain_id != target_chain:
                        continue
                    atom_name = line[12:16].strip()
                    if atom_name != 'CA':
                        continue
                    res_name = line[17:20].strip()
                    res_seq = line[22:26].strip()
                    i_code = line[26].strip()
                    imgt_pos = res_seq + i_code
                    aa = AA_MAP.get(res_name, 'X')
                    residues[imgt_pos] = aa
    except Exception as e:
        print(f"Error parsing {pdb_path}: {e}")
        return None
    return residues


def is_in_region(pos, start, end):
    num_part = "".join([c for c in pos if c.isdigit()])
    if not num_part:
        return False
    return start <= int(num_part) <= end


def sort_imgt_positions(positions):
    def key_func(pos):
        num_part = "".join([c for c in pos if c.isdigit()])
        char_part = "".join([c for c in pos if not c.isdigit()])
        return (int(num_part) if num_part else 0, char_part)
    return sorted(positions, key=key_func)


def extract_cdr(residues, start, end):
    positions = sort_imgt_positions([p for p in residues if is_in_region(p, start, end)])
    seq = "".join([residues[p] for p in positions])
    return len(positions), seq


def main():
    print("1. Loading metadata...")
    unique_pairs = load_unique_pairs(SUMMARY_FILE)
    print(f"   Found {len(unique_pairs)} unique (pdb, chain) pairs.")

    pdb_files = {os.path.basename(p)[:4].lower(): p
                 for p in glob.glob(os.path.join(IMGT_DIR, '*.pdb'))}
    print(f"   Found {len(pdb_files)} PDB files in '{IMGT_DIR}/'.")

    print("2. Calculating CDR lengths...")
    records = []
    skipped = 0
    for pdb_id, chain in unique_pairs:
        if pdb_id not in pdb_files:
            skipped += 1
            continue
        residues = parse_pdb_residues(pdb_files[pdb_id], chain)
        if residues is None:
            skipped += 1
            continue

        record = {'pdb_id': pdb_id, 'chain': chain}
        for cdr_name, (start, end) in CDR_REGIONS.items():
            length, seq = extract_cdr(residues, start, end)
            record[f'{cdr_name.lower()}_length'] = length
            record[f'{cdr_name.lower()}_sequence'] = seq if seq else None
        records.append(record)

    print(f"   Processed {len(records)} nanobodies, skipped {skipped}.")

    if not records:
        print("No valid data found.")
        return

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"   Saved results to {OUTPUT_CSV}")

    print("3. Generating visualization...")
    colors = {'CDR1': 'steelblue', 'CDR2': 'seagreen', 'CDR3': 'indianred'}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('CDR Length Distribution of Nanobody Structures (IMGT Numbering)',
                 fontsize=15, fontweight='bold', y=1.02)

    for ax, (cdr_name, (start, end)) in zip(axes, CDR_REGIONS.items()):
        col = f'{cdr_name.lower()}_length'
        lengths = df[col]
        n = len(lengths)
        mean_val = lengths.mean()
        median_val = lengths.median()
        std_val = lengths.std()

        length_counts = lengths.value_counts().sort_index()
        ax.bar(length_counts.index, length_counts.values,
               color=colors[cdr_name], edgecolor='white', linewidth=0.5)

        ax.axvline(mean_val, color='black', linestyle='--', linewidth=1.5,
                   label=f'Mean = {mean_val:.1f}')
        ax.axvline(median_val, color='dimgray', linestyle='-.', linewidth=1.5,
                   label=f'Median = {median_val:.1f}')

        ax.set_xlabel('Length (residues)', fontsize=12)
        ax.set_ylabel('Count (nanobodies)', fontsize=12)
        ax.set_title(f'{cdr_name} (IMGT pos {start}–{end})', fontsize=13)
        ax.legend(fontsize=10)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        stats_text = f'n = {n}\nmean = {mean_val:.1f}\nmedian = {median_val:.1f}\nstd = {std_val:.1f}'
        ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
                ha='right', va='top', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                          edgecolor='gray', alpha=0.8))

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_PLOT.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    print(f"   Saved plot to {OUTPUT_PLOT} / .pdf")
    print("Done.")


if __name__ == "__main__":
    main()
