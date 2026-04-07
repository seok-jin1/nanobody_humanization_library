#!/usr/bin/env python3
"""
IMGT 번호 체계를 PDB 순차 번호로 변환하여 Rosetta 돌연변이 파일 생성

Uses ANARCI to apply IMGT numbering to the nanobody sequence, then maps
each IMGT-numbered humanization mutation to the corresponding sequential
PDB position. Verifies wild-type residue identity at each mapped position.

Input:
    - (hardcoded) nanobody sequence: 125-residue anti-FAP VHH (CN106928368B)
    - (hardcoded) 13 IMGT-numbered humanization mutations
Output:
    - 01_mutations_final.txt : Rosetta mutation format file
      (each line: WT_residue position MUT_residue)
    - 01_imgt_pdb_mapping.csv : full IMGT ↔ PDB sequential position
      mapping table for all 125 residues

Usage:
    python 01_fix_imgt_numbering.py

Dependencies:
    anarci
"""

import csv
import os
import sys

try:
    import anarci
except ImportError:
    print("ERROR: anarci is required. Install with: pip install anarci")
    sys.exit(1)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

TARGET_SEQ = "QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS"

# 13 humanization mutations (IMGT numbering)
# Identified by aligning to IGHV3-66 via IMGT DomainGapAlign
# Hallmark residues excluded: F42, Q49, R50, A52
MUTATIONS_IMGT = [
    ("Q",   1, "E"),
    ("Q",   5, "V"),
    ("S",  12, "L"),
    ("A",  15, "P"),
    ("G",  40, "S"),
    ("V",  45, "A"),
    ("A",  54, "S"),
    ("I",  55, "V"),
    ("A",  83, "S"),
    ("K",  95, "R"),
    ("P",  96, "A"),
    ("M", 101, "V"),
    ("Q", 123, "L"),
]

OUTPUT_FILE = "01_mutations_final.txt"


IMGT_REGIONS = {
    "FR1": (1, 26), "CDR1": (27, 38), "FR2": (39, 55),
    "CDR2": (56, 65), "FR3": (66, 104), "CDR3": (105, 117), "FR4": (118, 128),
}


def get_region(pos):
    """Return IMGT region name for a given position number."""
    for name, (s, e) in IMGT_REGIONS.items():
        if s <= pos <= e:
            return name
    return "?"


def build_imgt_to_seq_map(sequence):
    """Use ANARCI to build IMGT position → sequential position mapping.

    Returns:
        imgt_to_seq: dict mapping IMGT main positions to (seq_pos, aa)
        full_mapping: list of (pdb_pos, aa, imgt_label, ins, region) for all residues
    """
    numbered, _, _ = anarci.anarci(
        [("nanobody", sequence)],
        scheme="imgt",
        allow={"H"},
        assign_germline=False,
    )

    if numbered is None or numbered[0] is None:
        print("ERROR: ANARCI failed to number the sequence.")
        sys.exit(1)

    domain, start, _ = numbered[0][0]

    imgt_to_seq = {}
    full_mapping = []
    seq_idx = start

    for (pos, ins), aa in domain:
        if aa == "-":
            continue
        seq_pos = seq_idx + 1  # 1-based PDB position
        ins_str = ins.strip()
        imgt_label = f"{pos}{ins_str}" if ins_str else str(pos)
        region = get_region(pos)

        full_mapping.append((seq_pos, aa, imgt_label, ins_str, region))

        if not ins_str:
            imgt_to_seq[pos] = (seq_pos, aa)
        seq_idx += 1

    return imgt_to_seq, full_mapping


def main():
    print(f"Target nanobody sequence ({len(TARGET_SEQ)} residues):")
    print(TARGET_SEQ)
    print()

    # Step 1: Build IMGT → sequential mapping via ANARCI
    print("Running ANARCI (IMGT scheme)...")
    imgt_to_seq, full_mapping = build_imgt_to_seq_map(TARGET_SEQ)
    print(f"  Mapped {len(imgt_to_seq)} IMGT positions to sequential positions.")
    print()

    # Step 2: Map and verify each mutation
    print("=" * 70)
    print("MUTATION MAPPING: IMGT → Sequential (PDB)")
    print("=" * 70)
    print(f"{'IMGT':>6}  {'WT':>3} → {'MUT':>3}  {'Seq Pos':>7}  {'Actual':>6}  {'Status'}")
    print("-" * 70)

    verified = []
    errors = []

    for wt, imgt_pos, mut in MUTATIONS_IMGT:
        if imgt_pos not in imgt_to_seq:
            print(f"{imgt_pos:>6}  {wt:>3} → {mut:>3}  {'N/A':>7}  {'N/A':>6}  ✗ IMGT position not found")
            errors.append((wt, imgt_pos, mut, "IMGT position not in ANARCI output"))
            continue

        seq_pos, actual_aa = imgt_to_seq[imgt_pos]

        if actual_aa == wt:
            print(f"{imgt_pos:>6}  {wt:>3} → {mut:>3}  {seq_pos:>7}  {actual_aa:>6}  ✓ MATCH")
            verified.append((wt, seq_pos, mut, imgt_pos))
        else:
            print(f"{imgt_pos:>6}  {wt:>3} → {mut:>3}  {seq_pos:>7}  {actual_aa:>6}  ✗ MISMATCH (expected {wt}, got {actual_aa})")
            errors.append((wt, imgt_pos, mut, f"Expected {wt} at seq {seq_pos}, found {actual_aa}"))

    print("-" * 70)
    print(f"Verified: {len(verified)} / {len(MUTATIONS_IMGT)}")

    if errors:
        print(f"\n⚠ ERRORS ({len(errors)}):")
        for wt, imgt_pos, mut, msg in errors:
            print(f"  IMGT {imgt_pos}: {wt}→{mut} — {msg}")
        sys.exit(1)

    # Step 3: Write Rosetta mutation file
    with open(OUTPUT_FILE, "w") as f:
        f.write(f"total {len(verified)}\n")
        for wt, seq_pos, mut, imgt_pos in verified:
            f.write("1\n")
            f.write(f"{wt} {seq_pos} {mut}\n")

    print(f"\n✓ Created: {OUTPUT_FILE}")
    print(f"  Contains {len(verified)} mutations in Rosetta format")

    # Step 4: Print mapping summary table
    print(f"\n{'IMGT':>6} {'Region':>6} {'WT':>3}→{'MUT':>3} {'SeqPos':>6}")
    print("-" * 35)
    for wt, seq_pos, mut, imgt_pos in verified:
        region = get_region(imgt_pos)
        print(f"{imgt_pos:>6} {region:>6} {wt:>3}→{mut:>3} {seq_pos:>6}")

    # Step 5: Save full IMGT ↔ PDB mapping table as CSV
    mutations_dict = {imgt: mut for _, imgt, mut in MUTATIONS_IMGT}
    mapping_csv = os.path.join(SCRIPT_DIR, "01_imgt_pdb_mapping.csv")
    with open(mapping_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_pos", "aa", "imgt_pos", "insertion", "region", "mutation"])
        for pdb_pos, aa, imgt_label, ins, region in full_mapping:
            # Check if this position has a humanization mutation
            imgt_num = int("".join(c for c in imgt_label if c.isdigit()))
            mut = ""
            if not ins and imgt_num in mutations_dict:
                mut = f"{aa}→{mutations_dict[imgt_num]}"
            writer.writerow([pdb_pos, aa, imgt_label, ins, region, mut])

    print(f"\n✓ Created: {mapping_csv}")
    print(f"  Full IMGT ↔ PDB mapping ({len(full_mapping)} residues)")


if __name__ == "__main__":
    main()
