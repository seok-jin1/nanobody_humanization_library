#!/usr/bin/env python3
"""
VH3_alleles_summary.csv 의 DNA sequence 컬럼을 아미노산 서열로 번역하여
aa_sequence 컬럼을 추가(또는 갱신)한다.

표준 유전암호(standard genetic code)를 사용하며,
인식 불가 코돈은 'X', 정지 코돈은 '*' 로 표기한다.
"""

import csv
import os

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, "VH3_alleles_summary.csv")


def translate(dna: str) -> str:
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3].upper()
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def main() -> None:
    rows = []
    with open(CSV_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            row['aa_sequence'] = translate(row['sequence'])
            rows.append(row)

    fieldnames = ['allele', 'gapped_length', 'ungapped_length', 'sequence', 'aa_sequence']
    with open(CSV_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done. {len(rows)} records updated: {CSV_PATH}")


if __name__ == "__main__":
    main()
