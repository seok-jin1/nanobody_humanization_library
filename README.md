# Humanized Nanobody Platform Construction

Computational pipeline for humanizing an anti-FAP nanobody (CN106928368B) using structure-guided framework engineering and CDR library design.

> Detailed research writeup with figures: [description.md](./description.md)

## Overview

```
Target: Fibroblast Activation Protein (FAP)
Source: CN106928368B anti-FAP nanobody (125 aa)
Closest human germline: IGHV3-66 (identified via IMGT DomainGapAlign)
Hallmark residues preserved: F42, Q49, R50, A52 (IMGT positions)
```

## Workflow

```
01_structure_prediction/     AlphaFold3 structure prediction
        │
        ▼
02_rosetta_ddg/              Rosetta Cartesian ΔΔG for 13 humanization mutations
        │                    → 7 confirmed + Set2 saturation (I55, P96)
        ▼
03_vh3_family_analysis/      Human VH3 germline AA proportion (IMGT/OGRDB)
        │                    → Residue-level conservation for borderline decisions
        ▼
04_sabdab_nano_analysis/     SAbDab-nano structural database analysis
        │                    → CDR length distribution & positional AA frequency
        ▼
05_gromacs_mmpbsa/           GROMACS MD + MM-PBSA validation (optional)
```

## Humanization Decision Summary

### Framework Mutations (13 candidates from IGHV3-66 alignment)

**Confirmed humanization mutations (8):**

| IMGT | Mutation | ΔΔG (REU) | Basis |
|------|----------|-----------|-------|
| 5 | Q5V | -0.09 | Neutral |
| 12 | S12L | -3.54 | Stabilizing |
| 40 | G40S | -1.38 | Stabilizing |
| 55 | I55L | -0.25 | Set2 saturation - only negative ΔΔG |
| 83 | A83S | -2.13 | Stabilizing |
| 95 | K95R | +0.23 | Neutral (conservative) |
| 101 | M101V | -2.27 | Stabilizing |
| 123 | Q123L | -0.59 | Stabilizing |

**Rejected (kept as original nanobody residue):**

| IMGT | Candidate | ΔΔG (REU) | Rejection Reason |
|------|-----------|-----------|------------------|
| 1 | Q1E | +2.09 | Destabilizing; Q present in 14.9% of VH3 alleles |
| 15 | A15P | +6.04 | Strongly destabilizing (backbone rigidity) |
| 45 | V45A | +2.76 | Destabilizing; A is 100% in VH3 but ΔΔG too high |
| 54 | A54S | +4.16 | Destabilizing; A present in 11.9% of VH3 alleles |
| 96 | P96A/V/D | +3.44/+3.68/+5.62 | All variants destabilizing |

### CDR Library Design (from SAbDab-nano analysis)

- **CDR1**: fixed length 8 aa
- **CDR2**: fixed length 8 aa
- **CDR3**: 3 length variants (14, 17, 21 aa) based on frequency peaks
- Position-specific AA diversity determined from SAbDab-nano residue frequency
- Library construction via trimer phosphoramidite primers

## Directory Structure

```
nanobody_humanization/
├── 01_structure_prediction/    # AF3 models (5), confidences, job request
├── 02_rosetta_ddg/             # ΔΔG calculations (Set1: 13 mutations, Set2: saturation)
├── 03_vh3_family_analysis/     # VH3 germline allele data from OGRDB (67 unique sequences)
├── 04_sabdab_nano_analysis/    # ~1,900 nanobody structures from SAbDab-nano
├── 05_gromacs_mmpbsa/          # MD simulation scripts (WT 100ns completed)
└── docs/                       # Summary, timeline, energy unit explanation
```

## Tools & Data Sources

| Tool | Purpose |
|------|---------|
| [AlphaFold3](https://alphafoldserver.com) | Nanobody 3D structure prediction |
| [Rosetta](https://www.rosettacommons.org) | Cartesian ΔΔG stability calculation |
| [IMGT DomainGapAlign](https://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi) | Germline alignment (→ IGHV3-66) |
| [OGRDB](https://ogrdb.airr-community.org) | Human VH3 germline allele sequences |
| [ANARCI](https://github.com/oxpig/ANARCI) | IMGT numbering assignment |
| [SAbDab-nano](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab) | Nanobody structure database |
| [GROMACS](https://www.gromacs.org) | Molecular dynamics simulation |
| [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) | MM-PBSA free energy calculation |

## References

1. Jumper, J., et al. (2024). AlphaFold3. *Nature*.
2. Kellogg, E.H., et al. (2011). Cartesian ddG. *Proteins*, 79(3), 830-838.
3. Lefranc, M.P., et al. (2015). IMGT numbering. *Dev. Comp. Immunol.*, 27(1), 55-77.
4. Dunbar, J., et al. (2014). SAbDab. *Nucleic Acids Res.*, 42(D1), D1140-D1146.
