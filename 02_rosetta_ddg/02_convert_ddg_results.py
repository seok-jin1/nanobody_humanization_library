#!/usr/bin/env python3
"""
Rosetta cartesian_ddg 원시 결과(.ddg)를 표준 ddG 예측 형식으로 변환

Parses the raw energy output from Rosetta's cartesian_ddg protocol,
calculates mean ddG values relative to the wild-type baseline with
propagated standard deviations, and writes standardized predictions.

Input:
    - <ddg_file> : Rosetta .ddg output file (e.g., mutations_final.ddg,
      mutations_set2.ddg) provided as CLI argument
Output:
    - 02_ddg_predictions.out : standardized ddG predictions
      (format: ddG: mutation_name mean_ddG std_dev)

Usage:
    python 02_convert_ddg_results.py mutations_final.ddg
"""

import sys
import re
from collections import defaultdict
import numpy as np

def parse_ddg_file(ddg_file):
    """Parse mutations_final.ddg and calculate ddG values."""
    energies = defaultdict(list)
    
    with open(ddg_file, 'r') as f:
        for line in f:
            if line.startswith('COMPLEX:'):
                # Parse line: COMPLEX:   Round#: MUT_NAME:  energy ...
                parts = line.split()
                round_info = parts[1]  # Round#:
                mut_name = parts[2].rstrip(':')  # MUT_NAME:
                energy = float(parts[3])  # energy value
                
                energies[mut_name].append(energy)
    
    return energies

def calculate_ddg(energies):
    """Calculate ddG values relative to WT."""
    if 'WT' not in energies:
        print("Error: No WT baseline found", file=sys.stderr)
        sys.exit(1)
    
    wt_energies = np.array(energies['WT'])
    wt_mean = np.mean(wt_energies)
    wt_std = np.std(wt_energies, ddof=1)
    
    results = {}
    results['WT'] = {
        'energy': wt_mean,
        'std': wt_std,
        'ddG': 0.0,
        'ddG_std': 0.0
    }
    
    for mut_name, mut_energies in energies.items():
        if mut_name == 'WT':
            continue
        
        mut_array = np.array(mut_energies)
        mut_mean = np.mean(mut_array)
        mut_std = np.std(mut_array, ddof=1)
        
        ddg = mut_mean - wt_mean
        # Propagate error: std(A-B) = sqrt(std(A)^2 + std(B)^2)
        ddg_std = np.sqrt(mut_std**2 + wt_std**2)
        
        results[mut_name] = {
            'energy': mut_mean,
            'std': mut_std,
            'ddG': ddg,
            'ddG_std': ddg_std
        }
    
    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python 02_convert_ddg_results.py mutations_final.ddg")
        sys.exit(1)
    
    ddg_file = sys.argv[1]
    
    # Parse raw energies
    energies = parse_ddg_file(ddg_file)
    
    # Calculate ddG values
    results = calculate_ddg(energies)
    
    # Output in expected format
    output_file = "02_ddg_predictions.out"
    with open(output_file, 'w') as f:
        # Header
        f.write("ddG: description total std\n")
        
        # WT baseline
        wt = results['WT']
        f.write(f"ddG: WT {wt['ddG']:.3f} {wt['ddG_std']:.3f}\n")
        
        # Mutations sorted by name
        for mut_name in sorted(results.keys()):
            if mut_name == 'WT':
                continue
            data = results[mut_name]
            f.write(f"ddG: {mut_name} {data['ddG']:.3f} {data['ddG_std']:.3f}\n")
    
    print(f"Converted {len(results)-1} mutations")
    print(f"Output written to: {output_file}")
    
    return results

if __name__ == "__main__":
    results = main()
