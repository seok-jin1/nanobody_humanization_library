# Humanized Nanobody Platform Construction

## 1. Overview

This study aims to establish a **humanized nanobody platform** for producing multi-specific antibodies targeting the tumor microenvironment (TME). Nanobodies (VHH), camelid-derived single-domain antibodies, offer several advantages over conventional IgG antibodies — small size (~15 kDa), high tissue penetration, and excellent stability — but can trigger immunogenicity in clinical applications. To address this, we established a humanization strategy that replaces framework (FR) regions with human germline sequences, and designed optimal sequences that maximize human similarity while preserving structural stability.

Building on this humanized FR scaffold, we construct nanobody libraries against novel targets by introducing diverse CDR sequences, and implement a platform that uses yeast display to select multi-specific antibodies binding to various targets within the tumor microenvironment (CAF markers, immune checkpoints, etc.).

---

## 2. Model Nanobody Selection and Experimental Validation

### 2.1 Anti-FAP Nanobody

An anti-Fibroblast Activation Protein (FAP) nanobody was selected as the model antibody (CN106928368B, 125 aa).

```
QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS
```

This nanobody is reliably produced in our lab via an *E. coli* expression system, showing high yield and stability after purification. It also exhibits **cross-reactivity** against human and mouse FAP, enabling efficacy evaluation in preclinical animal models.

### 2.2 Target Germline Selection for Humanization

Using IMGT DomainGapAlign (https://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi), **IGHV3-66** was identified as the human germline with the highest sequence similarity to the anti-FAP nanobody. The VH3 family is known to be structurally most similar to nanobodies (VHH) among human VH germlines, making it well-suited for CDR grafting.

Comparison with the IGHV3-66 germline revealed 13 differing residues in the FR regions, which were selected as humanization mutation candidates:

> Q1E, Q5V, S12L, A15P, G40S, V45A, A54S, I55V, A83S, K95R, P96A, M101V, Q123L

**VHH hallmark residues** (IMGT positions 42, 49, 50, 52), which play critical roles in nanobody solubility and stability, were excluded from humanization. By preserving the VHH-specific residues (F42, Q49, R50, A52) at these positions, the intrinsic structural features of the nanobody were maintained.

---

## 3. Structure-Based Stability Assessment

### 3.1 3D Structure Prediction

To assess the structural impact of humanization mutations, we first predicted the 3D structure of the nanobody. Using AlphaFold3, we obtained a high-confidence model of the anti-FAP nanobody (pTM = 0.88), and verified its suitability for downstream analysis by confirming disulfide bond formation (Cys22-Cys95), the absence of steric clashes, and the absence of disordered regions.

Per-residue pLDDT (predicted Local Distance Difference Test) analysis showed an overall average of 92.4, indicating very high confidence. FR regions predominantly showed very high confidence (>95), while only parts of the CDR3 loop ranged from 70-80 — a natural consequence of the loop's intrinsic flexibility. PAE (Predicted Aligned Error) analysis also showed uniformly low errors in the relative positioning of residues, confirming the reliability of the predicted structure.

![pLDDT Per-Residue](06_humanness_score/figures/04_plddt_per_residue.png)

![PAE Heatmap](06_humanness_score/figures/04_pae_heatmap.png)

### 3.2 Rosetta Cartesian ddG Calculation

Based on the predicted structure, we quantitatively assessed thermodynamic stability changes for each of the 13 candidate humanization mutations by calculating **Rosetta Cartesian ddG**.

Cartesian ddG computes the Rosetta energy difference between mutant and wild-type (ddG = E_mutant − E_WT) in Cartesian coordinate space, achieving higher accuracy than torsion-space methods (Kellogg et al., 2011; Park et al., 2016). The ref2015_cart score function was used, and 5 independent replicates were performed per mutation for statistical reliability. Key parameters are as follows:

- Score function: ref2015_cart (Cartesian-space optimization)
- Iterations: 5 (independent replicates)
- Backbone flexibility: ±1 residue (bbnbrs 1)
- Sidechain repacking radius: 9 Å (fa_max_dis 9.0)

Mutations were classified by ddG value as follows:

| Category | ddG range (REU) | Interpretation |
|------|:---:|------|
| Strongly destabilizing | > +2.0 | High risk to expression/stability |
| Moderately destabilizing | +1.0 to +2.0 | Moderate risk; consider reversion |
| Neutral | -0.5 to +0.5 | Minimal impact |
| Stabilizing | < -0.5 | Favorable for stability |

### 3.3 Set1 Results: Evaluation of 13 Individual Mutations

![Set1 ddG Bar Plot](02_rosetta_ddg/figures/04_ddg_barplot.png)

ddG analysis of the 13 mutations showed that 7 could be safely introduced as neutral or stabilizing, while 6 were classified as strongly destabilizing (ddG > +2.0 REU):

| Mutation | IMGT | ddG (REU) | Category |
|----------|:---:|:---:|------|
| Q5V | 5 | -0.09 | Neutral |
| S12L | 12 | -3.54 | Stabilizing |
| G40S | 40 | -1.38 | Stabilizing |
| A83S | 83 | -2.13 | Stabilizing |
| K95R | 95 | +0.23 | Neutral |
| M101V | 101 | -2.27 | Stabilizing |
| Q123L | 123 | -0.59 | Stabilizing |
| **Q1E** | **1** | **+2.09** | **Destabilizing** |
| **A15P** | **15** | **+6.04** | **Destabilizing** |
| **V45A** | **45** | **+2.76** | **Destabilizing** |
| **A54S** | **54** | **+4.16** | **Destabilizing** |
| **I55V** | **55** | **+2.52** | **Destabilizing** |
| **P96A** | **96** | **+3.44** | **Destabilizing** |

---

## 4. In-Depth Analysis of Destabilizing Mutations

For the 6 positions classified as destabilizing, two additional analyses were conducted: (1) amino acid frequency analysis within human VH3 germlines, and (2) exploration of alternative amino acids.

### 4.1 VH3 Germline Amino Acid Frequency Analysis

We collected 93 alleles (67 unique sequences) of the human IGHV3 family from OGRDB (Open Germline Receptor Database), applied IMGT numbering using ANARCI, and computed amino acid occurrence frequencies at each position.

![VH3 Proportion at Destabilizing Positions](02_rosetta_ddg/figures/06_vh3_proportion_destabilizing.png)

This analysis enabled the following decisions:

**IMGT 1 (Q→E):** E dominates in human VH3 at 83.6%, but the original nanobody residue Q is also present at 14.9%. Given the destabilizing ddG of +2.09 REU and Q's minority presence in VH3, we decided **not to mutate** this position.

**IMGT 54 (A→S):** S accounts for 73.1% in VH3, but the original residue A is also present at 11.9%. Given the strongly destabilizing ddG of +4.16 REU, we decided **not to mutate** this position.

**IMGT 15 (A→P) and IMGT 45 (V→A):** These residues are fully conserved in human germlines (P at 100% and A at 100%, respectively), but the severe structural destabilization (ddG of +6.04 and +2.76 REU) led us to decide **not to mutate** these positions.

**IMGT 55 (I→V) and IMGT 96 (P→A):** These two positions were selected for exploration of alternative amino acids.

### 4.2 Set2: Exploration of Alternative Amino Acids (Position Saturation)

For IMGT 55 and 96, we tested substitutions with various amino acids found in VH3 germlines.

**IMGT 55:** In VH3, 10 amino acids are evenly distributed (V 20.9%, G 19.4%, A 16.4%, etc.), making it worthwhile to explore alternatives. Among 8 variants tested (I55→G, A, S, Y, R, F, L, Q), **only I55L showed a negative ddG (−0.25 REU)**, enabling humanization while preserving stability.

**IMGT 96:** A (79.1%) and T (17.9%) predominate, so P96V and P96D were additionally tested, but both were destabilizing (P96V: +3.68, P96D: +5.62 REU). Consequently, we finalized the decision **not to mutate** this position.

---

## 5. Final Humanized Sequence Determination

![Decision Summary Table](02_rosetta_ddg/figures/05_decision_table.png)

Integrating the Rosetta ddG stability assessment with the VH3 germline frequency analysis, **8 of the 13 candidate mutations were finally adopted**:

| IMGT | Mutation | ddG (REU) | Decision | Rationale |
|:---:|---|:---:|:---:|------|
| 5 | Q→V | -0.09 | Accept | Neutral; 97.0% in VH3 |
| 12 | S→L | -3.54 | Accept | Stabilizing |
| 40 | G→S | -1.38 | Accept | Stabilizing |
| 55 | I→**L** | -0.25 | Accept | Only negative ddG in Set2 |
| 83 | A→S | -2.13 | Accept | Stabilizing |
| 95 | K→R | +0.23 | Accept | Neutral; conservative substitution |
| 101 | M→V | -2.27 | Accept | Stabilizing |
| 123 | Q→L | -0.59 | Accept | Stabilizing |

Five mutations were rejected due to structural destabilization risks: Q1E (+2.09), A15P (+6.04), V45A (+2.76), A54S (+4.16), and P96A (+3.44).

Final humanized sequence:

```
WT:        QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS
                *     *                       *              *                       *           *     *                           *
Humanized: QVQLVESGGGLVQAGGSLRLSCAASGYTVRSSYMSWFRQVPGKQREAVALITSGGTTYYADSVKGRFTISRDNSKNTLYLQMNSLRPEDTAVYYCAGRTGFIGGIWFRDRDYDYWGQGTLVTVSS
```

### 5.1 Quantitative Assessment of Humanization Effects

We quantified humanness before and after humanization using two metrics:

![Humanness Score Comparison](06_humanness_score/figures/01_humanness_comparison.png)

| Metric | WT | Humanized | Change |
|------|:---:|:---:|:---:|
| FR Identity (IGHV3-66) | 80.0% | **87.5%** | **+7.5%** |
| Position-Weighted Humanness (VH3 family) | 66.8% | **70.9%** | **+4.1%** |

FR Identity is the proportion of FR residues matching the target germline IGHV3-66; Position-Weighted Humanness is the mean per-position amino acid frequency across 67 unique VH3 family germlines. Both metrics improved meaningfully after humanization, while Rosetta ddG analysis confirmed that structural stability was preserved or improved.

AlphaFold3 prediction on the humanized sequence produced 5 models, all showing the same pTM (0.88) as WT, confirming no loss of structural prediction confidence due to humanization.

![AF3 pTM Comparison](06_humanness_score/figures/01_af3_ptm_comparison.png)

RMSD between WT and humanized structures was 0.169 Å (all-atom), 0.155 Å (backbone), and 0.158 Å (Cα) — all below 0.5 Å — quantitatively confirming that the 8 FR mutations have no substantive effect on the overall 3D structure of the nanobody.

![RMSD Comparison](06_humanness_score/figures/04_rmsd_comparison.png)

![Structure Overlay](06_humanness_score/figures/04_structure_overlay.png)

### 5.2 Multi-Mutant Stability Validation

Since individual ddG values evaluate only single substitutions, we performed **multi-mutant Cartesian ddG** calculations to validate the combined effect of introducing all 8 mutations simultaneously.

![Energy Decomposition](07_rosetta_humanized/figures/01_energy_decomposition.png)

With all 8 mutations applied simultaneously, ddG = **−7.710 ± 0.198 REU**, confirming that the humanized sequence is thermodynamically more stable than WT. Compared to the simple sum of individual ddG values (−10.02 REU), the epistatic effect is **+2.31 REU**, indicating weak inter-mutation interactions but still overall strongly stabilizing.

![Epistatic Effect](07_rosetta_humanized/figures/01_epistatic_effect.png)

Energy decomposition analysis showed fa_atr (van der Waals attractive) as the largest stabilization contributor at −9.30 REU, with fa_dun (rotamer energy, −3.40 REU) and p_aa_pp (Ramachandran probability, −1.46 REU) also contributing to stabilization.

Additionally, total-energy comparison using the **FastRelax** protocol showed Humanized (−447.89 ± 1.49 REU, n=12) with **−9.84 REU** lower total energy than WT (−438.05 ± 3.34 REU, n=5), independently confirming improved overall structural stability after humanization. The smaller standard deviation in Humanized versus WT also suggests a narrower, more stable energy landscape for the humanized sequence.

![Stability Comparison](07_rosetta_humanized/figures/03_stability_comparison.png)

---

## 6. Structural Database Analysis for CDR Library Design

To build a nanobody library by introducing diverse CDRs onto the humanized FR scaffold, we systematically analyzed ~1,900 nanobody structures registered in the **SAbDab-nano database**. The goal was to characterize CDR length distributions and per-position amino acid preferences observed in natural nanobodies, and thereby determine design parameters for a structurally valid CDR library.

### 6.1 CDR Length Distribution

![CDR Length Distribution](04_sabdab_nano_analysis/figures/05_cdr_length_distribution.png)

Length distributions were analyzed under IMGT numbering for CDR1 (positions 27-38), CDR2 (positions 56-65), and CDR3 (positions 105-117, including insertion codes).

- **CDR1:** predominantly fixed at **8 residues** (no gaps within IMGT positions 27-38)
- **CDR2:** predominantly fixed at **8 residues**
- **CDR3:** local maxima observed at **14, 17, and 21 residues**

Since CDR3 length diversity is a key determinant of nanobody antigen-binding specificity, we decided to design the library with **3 CDR3 length variants** (14, 17, and 21 residues). CDR1 and CDR2 were fixed at 8 aa each.

### 6.2 Per-Position Amino Acid Preferences

![Nanobody Residue Variability](04_sabdab_nano_analysis/figures/04_imgt_numeric_only_analysis.png)

Sequence diversity was quantified by calculating Shannon entropy at each IMGT position. CDR regions showed substantially higher entropy than FR regions, confirming that CDRs are the primary source of antigen-binding diversity. Per-position amino acid frequencies within CDRs were analyzed and used to determine the amino acid composition ratios at each position when designing trimer phosphoramidite primers.

### 6.3 Library Construction Strategy

Integrating the above analyses, we decided to synthesize the CDR library on the humanized FR scaffold using trimer phosphoramidite chemistry. Trimer phosphoramidite enables codon-level synthesis, excluding stop codons and allowing precise control of amino acid ratios at each position.

---

## 7. Conclusion

This study determined the optimal humanized sequence for the anti-FAP nanobody and established design parameters for the CDR library. Throughout humanization, we systematically leveraged computational approaches (AlphaFold3 structure prediction, Rosetta ddG stability assessment, and VH3 germline frequency analysis) to jointly optimize structural stability and human similarity, and we defined a library design strategy reflecting the CDR characteristics of natural nanobodies through SAbDab-nano database analysis.

This humanized nanobody platform is not limited to anti-FAP nanobodies; it represents infrastructure broadly applicable to the future development of nanobodies against diverse targets within the tumor microenvironment.
