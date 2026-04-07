#!/usr/bin/env python3
"""
Compare AlphaFold3 predicted structures with 1ZVH template-based model.
Calculate RMSD for backbone atoms (CA) and all heavy atoms.
"""

import subprocess
import os

# Paths
af3_dir = "/home/laugh/rosetta/alphafold3_prediction_FAP_Nb"
template_model = "relaxed/threaded_template_1ZVH_chainL_0001_0005.pdb"

print("="*80)
print("STRUCTURE COMPARISON: AlphaFold3 vs Template-Based Model")
print("="*80)
print()

# Convert CIF to PDB using Rosetta or pymol
# First, let's check if we have the tools
af3_models = []
for i in range(5):
    cif_file = f"{af3_dir}/fold_anti_fap_nb_model_{i}.cif"
    if os.path.exists(cif_file):
        af3_models.append((i, cif_file))
        print(f"Found AF3 model {i}: {cif_file}")

print(f"\nTemplate-based model: {template_model}")
print()

# Check confidence scores
print("="*80)
print("ALPHAFOLD3 CONFIDENCE SCORES")
print("="*80)

import json
for i in range(5):
    conf_file = f"{af3_dir}/fold_anti_fap_nb_summary_confidences_{i}.json"
    if os.path.exists(conf_file):
        with open(conf_file, 'r') as f:
            data = json.load(f)
            print(f"\nModel {i}:")
            for key, value in data.items():
                if isinstance(value, (int, float)):
                    print(f"  {key}: {value:.3f}")
                else:
                    print(f"  {key}: {value}")

print()
print("="*80)
print("STRUCTURE COMPARISON USING ROSETTA")
print("="*80)
print()

print("We can compare these structures using:")
print("1. Rosetta score_jd2 to calculate RMSD")
print("2. Visual inspection in PyMOL")
print("3. TM-align for structural alignment")
print()

# Create a script for Rosetta RMSD calculation
print("Creating RMSD calculation script...")

rmsd_xml = """
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref15" weights="ref2015.wts"/>
    </SCOREFXNS>

    <FILTERS>
        <Rmsd name="rmsd" threshold="999" confidence="0"/>
    </FILTERS>

    <PROTOCOLS>
        <Add filter="rmsd"/>
    </PROTOCOLS>

    <OUTPUT/>
</ROSETTASCRIPTS>
"""

with open('rmsd_calc.xml', 'w') as f:
    f.write(rmsd_xml)

print("✓ Created rmsd_calc.xml")
print()

# For now, let's just note what needs to be done
print("TO CALCULATE RMSD, RUN:")
print("-" * 80)
print()
print("# First, convert CIF to PDB (if needed)")
print("# AlphaFold3 CIF files can be read directly by PyMOL or converted")
print()
print("# Option 1: Using PyMOL")
print("pymol fold_anti_fap_nb_model_0.cif -c -d 'save model_0.pdb' -q")
print()
print("# Option 2: Using Rosetta score_jd2")
print("score_jd2.linuxgccrelease -in:file:s model_0.pdb -in:file:native template.pdb")
print()
print("# Option 3: Manual comparison in PyMOL")
print("pymol template.pdb model_0.cif")
print("PyMOL> align model_0, template")
print("PyMOL> rms_cur model_0 and name CA, template and name CA")
print()

print("="*80)
print("QUICK CHECK: Sequence identity")
print("="*80)
print()

# Check if sequences match
target_seq = "QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS"
print(f"Target sequence: {target_seq}")
print(f"Length: {len(target_seq)}")
print()
print("Both AlphaFold3 models and template-based model should have this sequence.")
print("The differences will be in the 3D coordinates (structure), not the sequence.")
print()

# Create a simple comparison script
comparison_script = """#!/bin/bash
# Structure comparison script

echo "Comparing AlphaFold3 models with template-based model..."
echo ""

# Copy AF3 models to working directory
cp /home/laugh/rosetta/alphafold3_prediction_FAP_Nb/fold_anti_fap_nb_model_*.cif .

# Launch PyMOL for visual comparison
echo "Launching PyMOL for visual comparison..."
echo ""
echo "In PyMOL, run:"
echo "  load relaxed/threaded_template_1ZVH_chainL_0001_0005.pdb, template"
echo "  load fold_anti_fap_nb_model_0.cif, af3_0"
echo "  load fold_anti_fap_nb_model_1.cif, af3_1"
echo "  align af3_0, template"
echo "  align af3_1, template"
echo "  rms_cur af3_0 and name CA, template and name CA"
echo ""
echo "This will show the structural overlay and RMSD."
"""

with open('compare_pymol.sh', 'w') as f:
    f.write(comparison_script)

os.chmod('compare_pymol.sh', 0o755)
print("✓ Created compare_pymol.sh")
print()

print("="*80)
print("EXPECTED DIFFERENCES")
print("="*80)
print("""
1. AlphaFold3:
   - AI-predicted structure (very accurate for well-studied proteins)
   - Multiple models showing prediction uncertainty
   - Confidence scores (pLDDT, pTM) indicate reliability

2. Template-based (1ZVH threading):
   - Based on homologous structure
   - Manually refined with Rosetta relaxation
   - Good for framework regions, less accurate for CDRs

3. Expected RMSD:
   - Framework regions: ~1-2 Å (should be similar)
   - CDR regions: ~3-5 Å (may differ significantly)
   - Overall C-alpha RMSD: ~2-3 Å expected

4. For ddG calculations:
   - Both structures should give similar stability predictions
   - Framework mutations (our focus) should be consistent
   - If RMSD < 3 Å, ddG results should be reliable
""")

print()
print("="*80)
print("NEXT STEPS")
print("="*80)
print("""
1. ✓ Compare AlphaFold3 models (check which has best confidence)
2. ⏳ Calculate RMSD between AF3 best model and template model
3. ⏳ If RMSD < 3 Å: Current ddG results are good
4. ⏳ If RMSD > 3 Å: Consider re-running ddG with AF3 structure
5. ⏳ Compare ddG results from both structures (if different)
""")

print()
print("Run './compare_pymol.sh' for visual comparison")
print("="*80)
