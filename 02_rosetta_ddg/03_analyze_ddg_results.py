#!/usr/bin/env python3
"""
ddG 결과를 분석하여 돌연변이 안정성 영향을 분류하고 권장사항 제시

Classifies each humanization mutation by its stability impact
(strongly destabilizing / moderately destabilizing / mildly destabilizing /
neutral / stabilizing) based on Cartesian ddG values, prints a detailed
console report with IMGT mapping and recommendations, and saves the
analysis summary.

Input:
    - <ddg_predictions_file> : standardized ddG predictions file
      (e.g., 02_ddg_predictions_final.out) provided as CLI argument
Output:
    - console : detailed stability analysis report
    - 03_ddg_analysis_report_<suffix>.txt : saved analysis summary
      suffix is derived from input filename (e.g., final, set2)

Usage:
    python 03_analyze_ddg_results.py 02_ddg_predictions_final.out  → 03_ddg_analysis_report_final.txt
    python 03_analyze_ddg_results.py 02_ddg_predictions_set2.out   → 03_ddg_analysis_report_set2.txt
"""

import re
import sys
from collections import defaultdict

def parse_ddg_results(ddg_file):
    """Parse ddG predictions output file."""
    results = {}
    current_mutation = None

    with open(ddg_file, 'r') as f:
        for line in f:
            # Look for mutation identifier
            if line.startswith('ddG:'):
                # Format: ddG: mutation_name ddG_value std_dev
                parts = line.split()
                if len(parts) >= 3:
                    mutation = parts[1]
                    try:
                        ddg_value = float(parts[2])
                        std_dev = float(parts[3]) if len(parts) > 3 else 0.0
                        results[mutation] = {
                            'ddG': ddg_value,
                            'std': std_dev
                        }
                    except ValueError:
                        continue

    return results

def classify_mutation(ddg, std):
    """Classify mutation based on ddG value."""
    if ddg > 2.0:
        return "STRONGLY DESTABILIZING", "❌"
    elif ddg > 1.0:
        return "MODERATELY DESTABILIZING", "⚠️"
    elif ddg > 0.5:
        return "MILDLY DESTABILIZING", "⚡"
    elif ddg > -0.5:
        return "NEUTRAL", "✓"
    else:
        return "STABILIZING", "✅"

def main():
    if len(sys.argv) < 2:
        print("Usage: python 03_analyze_ddg_results.py <02_ddg_predictions.out>")
        sys.exit(1)

    ddg_file = sys.argv[1]

    try:
        results = parse_ddg_results(ddg_file)
    except FileNotFoundError:
        print(f"Error: File {ddg_file} not found")
        sys.exit(1)

    if not results:
        print("No ddG results found in file")
        sys.exit(1)

    # IMGT mapping for reference (key = MUT_<seqpos><3-letter AA> from Rosetta output)
    imgt_mapping = {
        # Set1: 13 humanization mutations
        'MUT_1GLU':   (1,   'FR1'),  # Q1E
        'MUT_5VAL':   (5,   'FR1'),  # Q5V
        'MUT_11LEU':  (12,  'FR1'),  # S12L
        'MUT_14PRO':  (15,  'FR1'),  # A15P
        'MUT_35SER':  (40,  'FR2'),  # G40S
        'MUT_40ALA':  (45,  'FR2'),  # V45A
        'MUT_49SER':  (54,  'FR2'),  # A54S
        'MUT_50VAL':  (55,  'FR2'),  # I55V
        'MUT_74SER':  (83,  'FR3'),  # A83S
        'MUT_86ARG':  (95,  'FR3'),  # K95R
        'MUT_87ALA':  (96,  'FR3'),  # P96A
        'MUT_92VAL':  (101, 'FR3'),  # M101V
        'MUT_120LEU': (123, 'FR4'),  # Q123L
        # Set2: I55 saturation (seq pos 50)
        'MUT_50GLY':  (55,  'FR2'),  # I55G
        'MUT_50ALA':  (55,  'FR2'),  # I55A
        'MUT_50SER':  (55,  'FR2'),  # I55S
        'MUT_50TYR':  (55,  'FR2'),  # I55Y
        'MUT_50ARG':  (55,  'FR2'),  # I55R
        'MUT_50PHE':  (55,  'FR2'),  # I55F
        'MUT_50LEU':  (55,  'FR2'),  # I55L
        'MUT_50GLN':  (55,  'FR2'),  # I55Q
        # Set2: P96 saturation (seq pos 87)
        'MUT_87VAL':  (96,  'FR3'),  # P96V
        'MUT_87ASP':  (96,  'FR3'),  # P96D
    }

    print("="*80)
    print("NANOBODY HUMANIZATION - CARTESIAN ddG ANALYSIS")
    print("="*80)
    print()

    # Sort by ddG value (most destabilizing first)
    sorted_results = sorted(results.items(), key=lambda x: x[1]['ddG'], reverse=True)

    # Categorize mutations
    categories = defaultdict(list)
    for mut_name, data in sorted_results:
        if mut_name == 'WT':
            continue
        category, symbol = classify_mutation(data['ddG'], data['std'])
        categories[category].append((mut_name, data))

    # Print summary statistics
    print("SUMMARY STATISTICS")
    print("-" * 80)
    total_mutations = len([m for m in results if m != 'WT'])
    print(f"Total mutations analyzed: {total_mutations}")
    print()

    for category in ["STRONGLY DESTABILIZING", "MODERATELY DESTABILIZING",
                     "MILDLY DESTABILIZING", "NEUTRAL", "STABILIZING"]:
        count = len(categories[category])
        if count > 0:
            percentage = (count / total_mutations) * 100
            print(f"  {category:30s}: {count:2d} ({percentage:5.1f}%)")

    print()
    print("="*80)

    # Print detailed results by category
    for category in ["STRONGLY DESTABILIZING", "MODERATELY DESTABILIZING",
                     "MILDLY DESTABILIZING", "NEUTRAL", "STABILIZING"]:
        mutations_in_category = categories[category]
        if not mutations_in_category:
            continue

        _, symbol = classify_mutation(2.5 if "STRONGLY" in category else
                                     1.5 if "MODERATELY" in category else
                                     0.7 if "MILDLY" in category else
                                     0.0 if "NEUTRAL" in category else -1.0, 0)

        print(f"\n{symbol} {category}")
        print("-" * 80)

        for mut_name, data in mutations_in_category:
            ddg = data['ddG']
            std = data['std']

            # Get IMGT info
            imgt_info = ""
            if mut_name in imgt_mapping:
                imgt_pos, region = imgt_mapping[mut_name]
                imgt_info = f"IMGT {imgt_pos:3d} | {region:4s}"

            # Confidence indicator
            confidence = "HIGH" if std < 0.5 else "MED" if std < 1.0 else "LOW"

            print(f"  {mut_name:8s} | {imgt_info:15s} | ddG = {ddg:6.2f} ± {std:5.2f} REU | {confidence}")

    print()
    print("="*80)
    print("RECOMMENDATIONS")
    print("="*80)

    # High priority to revert
    high_risk = [m for m, d in sorted_results if d['ddG'] > 2.0 and m != 'WT']
    if high_risk:
        print("\n🚨 HIGH PRIORITY - Strongly recommend reverting:")
        for mut in high_risk:
            print(f"    - {mut}: ddG = {results[mut]['ddG']:.2f} REU")

    # Medium priority
    medium_risk = [m for m, d in sorted_results if 1.0 < d['ddG'] <= 2.0 and m != 'WT']
    if medium_risk:
        print("\n⚠️  MEDIUM PRIORITY - Consider reverting:")
        for mut in medium_risk:
            print(f"    - {mut}: ddG = {results[mut]['ddG']:.2f} REU")

    # Safe mutations
    safe = [m for m, d in sorted_results if d['ddG'] <= 0.5 and m != 'WT']
    if safe:
        print("\n✅ SAFE TO KEEP:")
        for mut in safe:
            print(f"    - {mut}: ddG = {results[mut]['ddG']:.2f} REU")

    print()
    print("="*80)
    print("IMPORTANT NOTES")
    print("="*80)
    print("""
1. P87A (IMGT 96): Proline removal may show artificially high ddG due to
   known Cartesian ddG artifacts. Interpret with caution.

2. A14P (IMGT 15): Proline introduction typically restricts backbone.
   High ddG expected and may be genuine stability concern.

3. G35S (IMGT 40): Loss of glycine flexibility. May impact local structure.

4. Mutations with high standard deviation (>1.0 REU) should be run with
   more iterations for better confidence.

5. Consider structural context: Surface mutations are more tolerant than
   buried core mutations.
    """)

    # Check WT reference
    if 'WT' in results:
        wt_ddg = results['WT']['ddG']
        print(f"WT reference ddG: {wt_ddg:.2f} REU")
        print("(Should be close to 0. Large deviation indicates calculation issues)")

    # Save report to file
    import io, contextlib
    buf = io.StringIO()
    # Re-run the printing to capture output
    with contextlib.redirect_stdout(buf):
        print("="*80)
        print("NANOBODY HUMANIZATION - CARTESIAN ddG ANALYSIS")
        print("="*80)
        print(f"\nTotal mutations analyzed: {total_mutations}\n")
        for category in ["STRONGLY DESTABILIZING", "MODERATELY DESTABILIZING",
                         "MILDLY DESTABILIZING", "NEUTRAL", "STABILIZING"]:
            muts = categories[category]
            if not muts:
                continue
            _, symbol = classify_mutation(2.5 if "STRONGLY" in category else
                                         1.5 if "MODERATELY" in category else
                                         0.7 if "MILDLY" in category else
                                         0.0 if "NEUTRAL" in category else -1.0, 0)
            print(f"\n{symbol} {category}")
            print("-" * 80)
            for mut_name, data in muts:
                imgt_info = ""
                if mut_name in imgt_mapping:
                    imgt_pos, region = imgt_mapping[mut_name]
                    imgt_info = f"IMGT {imgt_pos:3d} | {region:4s}"
                confidence = "HIGH" if data['std'] < 0.5 else "MED" if data['std'] < 1.0 else "LOW"
                print(f"  {mut_name:12s} | {imgt_info:15s} | ddG = {data['ddG']:6.2f} ± {data['std']:5.2f} REU | {confidence}")

    # Derive report filename from input: 02_ddg_predictions_final.out → 03_ddg_analysis_report_final.txt
    import os
    base = os.path.splitext(os.path.basename(ddg_file))[0]  # e.g. "02_ddg_predictions_final"
    suffix = base.replace("02_ddg_predictions_", "")  # e.g. "final" or "set2"
    report_file = f"03_ddg_analysis_report_{suffix}.txt"
    with open(report_file, 'w') as f:
        f.write(buf.getvalue())

    print()
    print("="*80)
    print(f"Analysis complete. Results saved to: {report_file}")
    print("="*80)

if __name__ == "__main__":
    main()
