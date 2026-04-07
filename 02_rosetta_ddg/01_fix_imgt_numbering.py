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

Usage:
    python 01_fix_imgt_numbering.py

Dependencies:
    anarci
"""

import sys

try:
    import anarci
except ImportError:
    print("ERROR: anarci is required. Install with: pip install anarci")
    sys.exit(1)

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


def build_imgt_to_seq_map(sequence):
    """Use ANARCI to build IMGT position → sequential position mapping."""
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

    # Build mapping: IMGT position (int) → (sequential_pos, amino_acid)
    # All non-gap residues increment sequential index (including insertions)
    imgt_to_seq = {}
    seq_idx = start  # ANARCI returns 0-based start index into the sequence

    for (pos, ins), aa in domain:
        if aa == "-":
            continue  # gap in IMGT numbering, no residue
        seq_pos = seq_idx + 1  # 1-based PDB position
        # Only store main positions (no insertion code) for framework mutations
        if not ins.strip():
            imgt_to_seq[pos] = (seq_pos, aa)
        seq_idx += 1  # always increment for non-gap residues

    return imgt_to_seq


def main():
    print(f"Target nanobody sequence ({len(TARGET_SEQ)} residues):")
    print(TARGET_SEQ)
    print()

    # Step 1: Build IMGT → sequential mapping via ANARCI
    print("Running ANARCI (IMGT scheme)...")
    imgt_to_seq = build_imgt_to_seq_map(TARGET_SEQ)
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
    regions = {range(1, 27): "FR1", range(27, 39): "CDR1", range(39, 56): "FR2",
               range(56, 66): "CDR2", range(66, 105): "FR3", range(105, 118): "CDR3",
               range(118, 129): "FR4"}
    for wt, seq_pos, mut, imgt_pos in verified:
        region = next((r for rng, r in regions.items() if imgt_pos in rng), "?")
        print(f"{imgt_pos:>6} {region:>6} {wt:>3}→{mut:>3} {seq_pos:>6}")


if __name__ == "__main__":
    main()
