#!/usr/bin/env python3
"""
FR 전체 saturation mutagenesis용 Rosetta mutation 파일 생성

Generates mutation files for exhaustive single-point saturation
mutagenesis across all framework (FR1-FR4) positions. Each FR residue
is mutated to all 19 other amino acids. CDR positions and Cys residues
(disulfide bond) are excluded.

Input:
    - 01_imgt_pdb_mapping.csv : IMGT ↔ PDB position mapping
Output:
    - saturation/mutations_fr_saturation.txt : single Rosetta mutation file (all mutations)
    - saturation/mutations_fr_saturation_split/ : individual mutation files for parallel execution
    - 12_fr_saturation_summary.csv : summary of all planned mutations

Usage:
    python 12_generate_saturation_mutations.py

Dependencies:
    None (standard library only)
"""

import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MAPPING_FILE = os.path.join(SCRIPT_DIR, "01_imgt_pdb_mapping.csv")
SAT_DIR = os.path.join(SCRIPT_DIR, "saturation")
SPLIT_DIR = os.path.join(SAT_DIR, "mutations_fr_saturation_split")

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

# Exclude Cys (disulfide bond C23-C104 in IMGT) from mutation targets
EXCLUDE_AA = {'C'}
# Only mutate FR positions
FR_REGIONS = {'FR1', 'FR2', 'FR3', 'FR4'}


def load_fr_positions(mapping_file):
    """Load FR positions from mapping CSV, excluding Cys and insertion positions."""
    positions = []
    with open(mapping_file) as f:
        for row in csv.DictReader(f):
            if row['region'] not in FR_REGIONS:
                continue
            if row['insertion']:  # skip insertion positions
                continue
            if row['aa'] in EXCLUDE_AA:
                continue
            positions.append({
                'pdb_pos': int(row['pdb_pos']),
                'aa': row['aa'],
                'imgt_pos': int(row['imgt_pos']),
                'region': row['region'],
            })
    return positions


def main():
    os.makedirs(SAT_DIR, exist_ok=True)
    os.makedirs(SPLIT_DIR, exist_ok=True)

    positions = load_fr_positions(MAPPING_FILE)
    print(f"FR positions (excluding Cys): {len(positions)}")

    # Generate all mutations
    all_mutations = []
    for pos in positions:
        wt = pos['aa']
        for mut_aa in AA_LIST:
            if mut_aa == wt:
                continue
            all_mutations.append({
                'pdb_pos': pos['pdb_pos'],
                'wt': wt,
                'mut': mut_aa,
                'imgt_pos': pos['imgt_pos'],
                'region': pos['region'],
            })

    print(f"Total mutations: {len(all_mutations)}")
    print(f"  FR1: {sum(1 for m in all_mutations if m['region'] == 'FR1')}")
    print(f"  FR2: {sum(1 for m in all_mutations if m['region'] == 'FR2')}")
    print(f"  FR3: {sum(1 for m in all_mutations if m['region'] == 'FR3')}")
    print(f"  FR4: {sum(1 for m in all_mutations if m['region'] == 'FR4')}")

    # Write single combined mutation file
    combined_file = os.path.join(SAT_DIR, "mutations_fr_saturation.txt")
    with open(combined_file, "w") as f:
        f.write(f"total {len(all_mutations)}\n")
        for m in all_mutations:
            f.write("1\n")
            f.write(f"{m['wt']} {m['pdb_pos']} {m['mut']}\n")
    print(f"\nSaved combined: {combined_file}")

    # Write individual mutation files (for parallel execution)
    # Group by position for efficient parallel runs
    positions_grouped = {}
    for m in all_mutations:
        key = m['pdb_pos']
        if key not in positions_grouped:
            positions_grouped[key] = []
        positions_grouped[key].append(m)

    for pdb_pos, muts in positions_grouped.items():
        region = muts[0]['region']
        imgt = muts[0]['imgt_pos']
        wt = muts[0]['wt']
        filename = f"mut_{region}_IMGT{imgt}_{wt}{pdb_pos}.txt"
        filepath = os.path.join(SPLIT_DIR, filename)
        with open(filepath, "w") as f:
            f.write(f"total {len(muts)}\n")
            for m in muts:
                f.write("1\n")
                f.write(f"{m['wt']} {m['pdb_pos']} {m['mut']}\n")

    print(f"Saved {len(positions_grouped)} split files: {SPLIT_DIR}/")

    # Summary CSV
    summary_csv = os.path.join(SCRIPT_DIR, "12_fr_saturation_summary.csv")
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "wt", "mut", "imgt_pos", "region", "mutation_label"])
        for m in all_mutations:
            label = f"{m['wt']}{m['imgt_pos']}{m['mut']}"
            writer.writerow([m['pdb_pos'], m['wt'], m['mut'], m['imgt_pos'], m['region'], label])
    print(f"Saved summary: {summary_csv}")

    # Estimate computation time
    n_calc = (len(all_mutations) + len(positions_grouped)) * 5  # mutations + WT per file × iterations
    time_per_calc_min = 18 / 14  # ~1.3 min per calculation (based on Set1: 14 structs × 5 iter = 4hr)
    total_hours = n_calc * time_per_calc_min / 60
    print(f"\nEstimated computation:")
    print(f"  Calculations: ~{n_calc}")
    print(f"  Single core: ~{total_hours:.0f} hours")
    print(f"  6 cores: ~{total_hours / 6:.0f} hours (~{total_hours / 6 / 24:.1f} days)")


if __name__ == "__main__":
    main()
