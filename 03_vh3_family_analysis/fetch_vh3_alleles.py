#!/usr/bin/env python3
"""
Fetch human VH3 family (IGHV3-*) germline allele sequences from OGRDB REST API.
Saves gapped FASTA, ungapped FASTA, and summary CSV to the same directory.
"""

import csv
import os
import sys
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

OGRDB_URL = "https://ogrdb.airr-community.org/api_v2/germline/set/9606.IGH_VDJ/latest/gapped"

OUT_GAPPED   = os.path.join(SCRIPT_DIR, "VH3_alleles_gapped.fasta")
OUT_UNGAPPED = os.path.join(SCRIPT_DIR, "VH3_alleles_ungapped.fasta")
OUT_CSV      = os.path.join(SCRIPT_DIR, "VH3_alleles_summary.csv")


def fetch_url(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "fetch_vh3_alleles/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Return list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def filter_vh3(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return [(h, s) for h, s in records if h.startswith("IGHV3-")]


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def write_csv(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["allele", "gapped_length", "ungapped_length", "sequence"])
        for header, seq in records:
            ungapped = seq.replace(".", "").replace("-", "")
            writer.writerow([header, len(seq), len(ungapped), ungapped])


def main() -> None:
    print(f"Fetching from: {OGRDB_URL}")
    try:
        raw = fetch_url(OGRDB_URL)
    except urllib.error.URLError as e:
        print(f"OGRDB fetch failed: {e}", file=sys.stderr)
        sys.exit(1)

    all_records = parse_fasta(raw)
    print(f"Total IGHV records fetched: {len(all_records)}")

    vh3_records = filter_vh3(all_records)
    print(f"VH3 family records (IGHV3-*): {len(vh3_records)}")

    if not vh3_records:
        print("No VH3 records found. Check API response.", file=sys.stderr)
        sys.exit(1)

    # Gapped FASTA
    write_fasta(OUT_GAPPED, vh3_records)
    print(f"Saved gapped FASTA  : {OUT_GAPPED}")

    # Ungapped FASTA
    ungapped_records = [(h, s.replace(".", "").replace("-", "")) for h, s in vh3_records]
    write_fasta(OUT_UNGAPPED, ungapped_records)
    print(f"Saved ungapped FASTA: {OUT_UNGAPPED}")

    # CSV summary
    write_csv(OUT_CSV, vh3_records)
    print(f"Saved CSV summary   : {OUT_CSV}")

    print("\nDone.")


if __name__ == "__main__":
    main()
