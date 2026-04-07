#!/usr/bin/env python3
"""
VH3_alleles_summary.csv 에서 unique 아미노산 서열을 추출하고
ANARCI 를 통해 IMGT numbering 을 적용한다.

출력:
  VH3_unique_imgt.csv     — unique 서열별 IMGT position 컬럼 CSV
  VH3_unique.fasta        — unique 서열 FASTA
"""

import csv
import os
from collections import defaultdict

import anarci

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
IN_CSV       = os.path.join(SCRIPT_DIR, "VH3_alleles_summary.csv")
OUT_CSV      = os.path.join(SCRIPT_DIR, "VH3_unique_imgt.csv")
OUT_FASTA    = os.path.join(SCRIPT_DIR, "VH3_unique.fasta")


# ── 1. unique 서열 추출 ────────────────────────────────────────────────────────

def load_unique(csv_path: str) -> list[dict]:
    """
    aa_sequence 기준으로 중복 제거.
    같은 서열을 가진 allele 이름은 모두 보존 (all_alleles 필드).
    대표 이름(representative)은 처음 등장한 allele.
    """
    seen: dict[str, list[str]] = defaultdict(list)  # seq -> [allele, ...]
    order: list[str] = []  # 등장 순서 보존

    with open(csv_path) as f:
        for row in csv.DictReader(f):
            seq = row['aa_sequence']
            allele = row['allele']
            if seq not in seen:
                order.append(seq)
            seen[seq].append(allele)

    records = []
    for seq in order:
        alleles = seen[seq]
        records.append({
            'representative': alleles[0],
            'all_alleles': ';'.join(alleles),
            'aa_sequence': seq,
        })
    return records


# ── 2. IMGT numbering ─────────────────────────────────────────────────────────

def imgt_label(pos: int, ins: str) -> str:
    """(112, 'A') -> '112A',  (5, ' ') -> '5'"""
    return f"{pos}{ins.strip()}"


def run_anarci(records: list[dict]) -> list[dict]:
    """
    ANARCI 로 IMGT numbering 수행.
    각 record 에 'numbering' dict (position_label -> aa) 추가.
    numbering 실패 시 'numbering' = None.
    """
    seqs = [(r['representative'], r['aa_sequence']) for r in records]
    numbered, details, _ = anarci.anarci(
        seqs,
        scheme='imgt',
        allow={'H'},
        assign_germline=False,
    )

    for rec, num_result, detail in zip(records, numbered, details):
        if num_result is None or detail is None:
            rec['numbering'] = None
            rec['chain_type'] = None
            continue
        # num_result: list of (domain_numbering, start, end)
        # 첫 번째 도메인만 사용
        domain, start, end = num_result[0]
        rec['numbering'] = {
            imgt_label(pos, ins): aa
            for (pos, ins), aa in domain
        }
        rec['chain_type'] = detail[0]['chain_type'] if detail else None

    return records


# ── 3. 모든 IMGT position 컬럼 수집 ───────────────────────────────────────────

def collect_positions(records: list[dict]) -> list[str]:
    """
    전체 record 에 걸쳐 등장하는 IMGT position 을 순서대로 반환.
    정렬 기준: (숫자, 삽입 문자) — IMGT 순서와 동일.
    """
    pos_set: set[str] = set()
    for r in records:
        if r['numbering']:
            pos_set.update(r['numbering'].keys())

    def sort_key(label: str):
        digits = ''.join(c for c in label if c.isdigit())
        letters = ''.join(c for c in label if c.isalpha())
        return (int(digits), letters)

    return sorted(pos_set, key=sort_key)


# ── 4. 저장 ───────────────────────────────────────────────────────────────────

def write_csv(path: str, records: list[dict], positions: list[str]) -> None:
    meta_cols = ['representative', 'all_alleles', 'chain_type', 'aa_sequence']
    fieldnames = meta_cols + positions

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for r in records:
            row = {
                'representative': r['representative'],
                'all_alleles':    r['all_alleles'],
                'chain_type':     r.get('chain_type', ''),
                'aa_sequence':    r['aa_sequence'],
            }
            if r['numbering']:
                for pos in positions:
                    aa = r['numbering'].get(pos, '-')
                    row[pos] = aa
            else:
                for pos in positions:
                    row[pos] = ''
            writer.writerow(row)


def write_fasta(path: str, records: list[dict]) -> None:
    with open(path, 'w') as f:
        for r in records:
            f.write(f">{r['representative']}\n")
            seq = r['aa_sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("Loading CSV and extracting unique sequences...")
    records = load_unique(IN_CSV)
    print(f"  Total alleles : {sum(len(r['all_alleles'].split(';')) for r in records)}")
    print(f"  Unique seqs   : {len(records)}")

    print("Running ANARCI (IMGT scheme)...")
    records = run_anarci(records)
    failed = [r['representative'] for r in records if r['numbering'] is None]
    if failed:
        print(f"  WARNING: ANARCI failed for {len(failed)} sequences: {failed}")
    else:
        print(f"  All {len(records)} sequences numbered successfully.")

    positions = collect_positions(records)
    print(f"  IMGT positions : {positions[0]} ~ {positions[-1]}  ({len(positions)} columns)")

    write_csv(OUT_CSV, records, positions)
    print(f"Saved: {OUT_CSV}")

    write_fasta(OUT_FASTA, records)
    print(f"Saved: {OUT_FASTA}")


if __name__ == '__main__':
    main()
