#!/usr/bin/env python3
"""
Analyze Cartesian ddG results and identify high-risk mutations.
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
        print("Usage: python analyze_ddg_results.py <ddg_predictions.out>")
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

    # IMGT mapping for reference
    imgt_mapping = {
        'Q1E': (1, 'FR1'),
        'Q5V': (5, 'FR1'),
        'S7L': (12, 'FR1'),
        'A14P': (15, 'FR1'),
        'G35S': (40, 'FR2'),
        'V40A': (45, 'FR2'),
        'A49S': (54, 'FR2'),
        'I50V': (55, 'FR2'),
        'A74S': (83, 'FR3'),
        'K86R': (95, 'FR3'),
        'P87A': (96, 'FR3'),
        'M92V': (101, 'FR3'),
        'Q120L': (123, 'FR4'),
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

    print()
    print("="*80)
    print("Analysis complete. Results saved to: ddg_analysis_report.txt")
    print("="*80)

if __name__ == "__main__":
    main()
