#!/usr/bin/env python3
"""
IMGT position별 아미노산 빈도 계산

Calculates amino acid proportions at each IMGT position.
Gaps ('-') are excluded from proportion calculations.

Input:
    - 03_VH3_unique_imgt.csv : IMGT 넘버링이 적용된 unique 서열 CSV (03_imgt_numbering.py 출력)
Output:
    - 04_VH3_imgt_aa_proportion.csv      : wide format (rows=position, cols=amino acid)
    - 04_VH3_imgt_aa_proportion_long.csv : long format (position, amino_acid, count, proportion)

Usage:
    python 04_aa_proportion.py
"""

import csv
import os
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_CSV     = os.path.join(SCRIPT_DIR, "03_VH3_unique_imgt.csv")
OUT_WIDE   = os.path.join(SCRIPT_DIR, "04_VH3_imgt_aa_proportion.csv")
OUT_LONG   = os.path.join(SCRIPT_DIR, "04_VH3_imgt_aa_proportion_long.csv")

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY") + ["X", "*"]


def sort_key(label: str):
    digits = ''.join(c for c in label if c.isdigit())
    letters = ''.join(c for c in label if c.isalpha())
    return (int(digits), letters)


def main():
    with open(IN_CSV) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    pos_cols = sorted(
        [c for c in fieldnames if c[0].isdigit()],
        key=sort_key
    )
    n_seq = len(rows)
    print(f"Sequences : {n_seq}")
    print(f"Positions : {len(pos_cols)}  ({pos_cols[0]} ~ {pos_cols[-1]})")

    # position별 집계 (gap '-' 제외)
    stats: list[dict] = []
    all_aas: set[str] = set()

    for pos in pos_cols:
        residues = [r[pos] for r in rows if r[pos] not in ('-', '')]
        counter = Counter(residues)
        n_obs = sum(counter.values())       # gap 제외 관측 수
        all_aas.update(counter.keys())
        stats.append({
            'position': pos,
            'n_obs': n_obs,
            'gap_count': n_seq - n_obs,
            'gap_fraction': round((n_seq - n_obs) / n_seq, 4),
            'counter': counter,
        })

    # 등장한 아미노산만 컬럼으로 (AA_ORDER 기준 정렬)
    present_aas = [aa for aa in AA_ORDER if aa in all_aas]
    extra_aas   = sorted(all_aas - set(AA_ORDER))
    aa_cols = present_aas + extra_aas

    # ── wide format ───────────────────────────────────────────────────────────
    wide_header = ['position', 'n_obs', 'gap_count', 'gap_fraction'] + aa_cols
    with open(OUT_WIDE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=wide_header)
        writer.writeheader()
        for s in stats:
            row = {
                'position':     s['position'],
                'n_obs':        s['n_obs'],
                'gap_count':    s['gap_count'],
                'gap_fraction': s['gap_fraction'],
            }
            n = s['n_obs'] if s['n_obs'] > 0 else 1
            for aa in aa_cols:
                cnt = s['counter'].get(aa, 0)
                row[aa] = round(cnt / n, 4)
            writer.writerow(row)
    print(f"Saved wide : {OUT_WIDE}")

    # ── long format ───────────────────────────────────────────────────────────
    with open(OUT_LONG, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['position', 'amino_acid', 'count', 'proportion'])
        for s in stats:
            n = s['n_obs'] if s['n_obs'] > 0 else 1
            for aa, cnt in sorted(s['counter'].items(),
                                  key=lambda x: -x[1]):
                writer.writerow([
                    s['position'],
                    aa,
                    cnt,
                    round(cnt / n, 4),
                ])
    print(f"Saved long : {OUT_LONG}")

    # ── 콘솔 요약: 가장 다양한 position top 10 ────────────────────────────────
    diversity = [(s['position'], len(s['counter'])) for s in stats]
    diversity.sort(key=lambda x: -x[1])
    print("\nTop 10 most variable IMGT positions:")
    print(f"  {'pos':>6}  {'#aa':>4}  top residues")
    for pos, n_aa in diversity[:10]:
        s = next(x for x in stats if x['position'] == pos)
        top = ', '.join(
            f"{aa}:{s['counter'][aa]}"
            for aa, _ in s['counter'].most_common(3)
        )
        print(f"  {pos:>6}  {n_aa:>4}  {top}")


if __name__ == '__main__':
    main()
