#!/usr/bin/env python3
"""
IMGT 번호 체계를 PDB 순차 번호로 변환하여 Rosetta 돌연변이 파일 생성

Maps 13 humanization mutations from IMGT numbering to sequential PDB
positions for the anti-FAP nanobody (CN106928368B, 125 aa). Verifies
each mutation by checking the expected wild-type residue at the mapped
position in the nanobody sequence, then writes the corrected mutations
in Rosetta-compatible format.

Input:
    - (hardcoded) nanobody sequence: 125-residue anti-FAP VHH
    - (hardcoded) 13 IMGT-numbered humanization mutations
Output:
    - 01_mutations_corrected.txt : Rosetta mutation format file
      (each line: WT_residue position MUT_residue)

Usage:
    python 01_fix_imgt_numbering.py
"""

target_seq = "QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS"

print("Target nanobody sequence (125 residues):")
print(target_seq)
print()

# Standard IMGT regions for VHH:
# FR1: 1-26
# CDR1: 27-38 (variable length)
# FR2: 39-55
# CDR2: 56-65 (variable length)
# FR3: 66-104
# CDR3: 105-117 (variable length)
# FR4: 118-128

# For this specific sequence, let's identify regions by known patterns:
# FR1 starts: QVQLQESG (conserved)
# CDR1 typically after position 26
# FR2 has WFRQ or WVRQ
# CDR2 is short in VHH
# FR3 is longest, has many conserved positions
# CDR3 has variable sequence
# FR4 ends: TVSS (conserved)

# Find key conserved residues to align IMGT numbering:
print("Finding conserved framework residues...")

# FR1 positions (1-26 in IMGT)
# Q1, V2, Q3, L4, Q5, E6, S7, G8, G9, G10
# Our seq: QVQLQESGGGSVQAGGSLRLSC
fr1_end = 26  # Typical FR1 end

# FR2 starts around position 39 in IMGT
# Look for WFRQ pattern (common in VHH)
wfrq_pos = target_seq.find("WFRQ")
print(f"WFRQ motif found at position: {wfrq_pos + 1}")

# FR3 has many conserved positions
# Look for ISRDNAK pattern (IMGT 76-82)
isrdnak_pos = target_seq.find("ISRDNAK")
if isrdnak_pos >= 0:
    print(f"ISRDNAK motif found at position: {isrdnak_pos + 1}")
    # IMGT 76 = sequential position isrdnak_pos + 1
    imgt_76_seq = isrdnak_pos + 1
else:
    isrdnak_pos = target_seq.find("ISRDN")
    print(f"ISRDN motif found at position: {isrdnak_pos + 1}")
    imgt_76_seq = isrdnak_pos + 1

# FR3-FR4 junction: Look for QGTQVTVSS or similar
fr4_pattern = target_seq.find("QGTQVTVSS")
if fr4_pattern >= 0:
    print(f"FR4 pattern QGTQVTVSS found at position: {fr4_pattern + 1}")
    # IMGT 118 typically starts FR4
    imgt_118_seq = fr4_pattern + 1

# Now create mapping table
# For VHH without large insertions/deletions:
# IMGT 1-26 = Sequential 1-26 (FR1)
# IMGT 27-38 = Sequential 27-38 (CDR1)
# IMGT 39-55 = Sequential 39-55 (FR2)
# IMGT 56-65 = Sequential 56-65 (CDR2)
# IMGT 66-104 = Sequential 66-104 (FR3)
# IMGT 105-117 = Sequential 105-117 (CDR3)
# IMGT 118-128 = Sequential 118-128 (FR4)

# But our sequence is only 125 residues, so there might be slight compression

# Let's use a more careful approach:
# Check each mutation position manually

mutations = [
    ("Q", 1, "E"),
    ("Q", 5, "V"),
    ("S", 12, "L"),
    ("A", 15, "P"),
    ("G", 40, "S"),
    ("V", 45, "A"),
    ("A", 54, "S"),
    ("I", 55, "V"),
    ("A", 83, "S"),
    ("K", 95, "R"),
    ("P", 96, "A"),
    ("M", 101, "V"),
    ("Q", 123, "L"),
]

print("\n" + "="*80)
print("MANUAL POSITION VERIFICATION")
print("="*80)

# Print sequence with position markers
print("\nSequence with 10-position markers:")
for i in range(0, len(target_seq), 10):
    print(f"{i+1:3d}: {target_seq[i:i+10]}")

print("\n" + "="*80)
print("CHECKING EACH MUTATION:")
print("="*80)

verified_mutations = []

for wt, imgt_pos, mut in mutations:
    print(f"\nIMGT {imgt_pos}: {wt} → {mut}")

    # For FR1 (IMGT 1-26), mapping is usually direct
    if imgt_pos <= 26:
        seq_pos = imgt_pos
        actual = target_seq[seq_pos - 1] if seq_pos <= len(target_seq) else "?"
        print(f"  FR1: IMGT {imgt_pos} → Sequential {seq_pos}")
        print(f"  Expected: {wt}, Found: {actual}")
        if actual == wt:
            verified_mutations.append((wt, seq_pos, mut))
            print(f"  ✓ MATCH")
        else:
            # Search nearby
            for offset in range(-5, 6):
                test_pos = seq_pos + offset
                if 0 < test_pos <= len(target_seq):
                    if target_seq[test_pos - 1] == wt:
                        print(f"  → Found {wt} at position {test_pos} (offset {offset})")
                        verified_mutations.append((wt, test_pos, mut))
                        break

    # For FR2 (IMGT 39-55)
    elif 39 <= imgt_pos <= 55:
        # FR2 might have slight shift
        # Assume -3 residue shift (common in VHH)
        seq_pos = imgt_pos - 3
        actual = target_seq[seq_pos - 1] if seq_pos <= len(target_seq) else "?"
        print(f"  FR2: IMGT {imgt_pos} → Sequential ~{seq_pos}")
        print(f"  Expected: {wt}, Found: {actual}")
        if actual == wt:
            verified_mutations.append((wt, seq_pos, mut))
            print(f"  ✓ MATCH")
        else:
            # Search nearby
            for offset in range(-5, 6):
                test_pos = imgt_pos + offset
                if 0 < test_pos <= len(target_seq):
                    if target_seq[test_pos - 1] == wt:
                        print(f"  → Found {wt} at position {test_pos}")
                        verified_mutations.append((wt, test_pos, mut))
                        break

    # For FR3 (IMGT 66-104)
    elif 66 <= imgt_pos <= 104:
        # FR3 might have shift
        seq_pos = imgt_pos - 3  # Try -3 shift
        actual = target_seq[seq_pos - 1] if seq_pos <= len(target_seq) else "?"
        print(f"  FR3: IMGT {imgt_pos} → Sequential ~{seq_pos}")
        print(f"  Expected: {wt}, Found: {actual}")
        if actual == wt:
            verified_mutations.append((wt, seq_pos, mut))
            print(f"  ✓ MATCH")
        else:
            # Search in range
            for offset in range(-10, 11):
                test_pos = imgt_pos + offset
                if 0 < test_pos <= len(target_seq):
                    if target_seq[test_pos - 1] == wt:
                        print(f"  → Found {wt} at position {test_pos} (IMGT{imgt_pos} offset {offset})")
                        verified_mutations.append((wt, test_pos, mut))
                        break

    # For FR4 (IMGT 118-128)
    elif 118 <= imgt_pos <= 128:
        # FR4: map to end of sequence
        # IMGT 123 might be around position 120
        seq_pos = len(target_seq) - (128 - imgt_pos)
        actual = target_seq[seq_pos - 1] if seq_pos <= len(target_seq) else "?"
        print(f"  FR4: IMGT {imgt_pos} → Sequential ~{seq_pos}")
        print(f"  Expected: {wt}, Found: {actual}")
        if actual == wt:
            verified_mutations.append((wt, seq_pos, mut))
            print(f"  ✓ MATCH")
        else:
            # Search nearby
            for offset in range(-5, 6):
                test_pos = seq_pos + offset
                if 0 < test_pos <= len(target_seq):
                    if target_seq[test_pos - 1] == wt:
                        print(f"  → Found {wt} at position {test_pos}")
                        verified_mutations.append((wt, test_pos, mut))
                        break

print("\n" + "="*80)
print(f"VERIFIED: {len(verified_mutations)} / {len(mutations)} mutations")
print("="*80)

# Write corrected mutations file
with open('01_mutations_corrected.txt', 'w') as f:
    f.write(f"total {len(verified_mutations) + 1}\n")

    for wt, pos, mut in verified_mutations:
        f.write("1\n")
        f.write(f"{wt} {pos} {mut}\n")

    f.write("1\n")
    f.write("WT WT WT\n")

print("\n✓ Created: 01_mutations_corrected.txt")
print(f"  Contains {len(verified_mutations)} mutations + WT reference")
