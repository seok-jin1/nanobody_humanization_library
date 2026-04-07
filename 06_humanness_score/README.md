# Humanness Score & Structural Comparison

인간화 전후 nanobody의 humanness를 정량적으로 평가하고,
GROMACS MD 시뮬레이션으로 구조적 안정성을 비교한다.

---

## 배경

Nanobody(VHH)는 낙타과 유래 단일 도메인 항체로, 인체 투여 시 면역원성 위험이 있다.
인간화(humanization)는 CDR을 유지하면서 framework(FR) 잔기를 human germline에 가깝게 변환하여
면역원성을 낮추는 전략이다.

VHH에 특화된 humanness 지표는 아직 확립되지 않았으며(T20, OASis, Hu-mAb 등은 모두
conventional VH/VL 기준), 대부분의 nanobody 인간화 논문에서는 **가장 가까운 human germline과의
FR identity**를 humanness 지표로 사용한다 (Vincke et al., 2009; Rossotti et al., 2022).

---

## 지표

### 1. FR Identity to IGHV3-66

FR 잔기 중 target germline (IGHV3-66\*01)과 동일한 잔기의 비율.

```
FR Identity = (IGHV3-66과 일치하는 FR 잔기 수) / (비교 가능한 FR 잔기 수)
```

- 가장 직관적이고 논문에서 가장 널리 사용되는 지표
- CDR 및 VHH hallmark 잔기(IMGT 42, 49, 50, 52)는 의도적으로 변경하지 않았으므로 FR만 비교

### 2. Position-Weighted Humanness (PWH)

각 FR position에서 해당 AA가 VH3 family 전체(67 unique germline)에서 차지하는 비율의 평균.

```
PWH = mean(VH3_proportion(AA_i) for each FR position i)
```

- FR Identity가 단일 germline(IGHV3-66) 기준인 반면, PWH는 VH3 family 전체를 반영
- 특정 position에서 germline과 불일치하더라도 VH3 내에서 흔한 AA이면 높은 점수

---

## 결과 요약

| 지표 | WT (original) | Humanized | 변화 |
|------|:---:|:---:|:---:|
| FR Identity (IGHV3-66) | 80.0% (64/80) | 87.5% (70/80) | **+7.5%** |
| Position-Weighted Humanness (VH3) | 66.8% | 70.9% | **+4.1%** |

8개 mutation 적용: Q5V, S12L, G40S, I55L, A83S, K95R, M101V, Q123L

---

## 파일 목록

### 스크립트

| 파일 | 설명 |
|------|------|
| `01_humanness_score.py` | humanness score 계산 + 시각화 |
| `02_run_md_humanized.sh` | Humanized nanobody 100ns MD 실행 (GROMACS) |
| `03_extract_xvg.sh` | WT/Humanized trajectory에서 RMSD, RMSF, Rg 추출 |
| `03_compare_md.py` | MD 결과 비교 분석 + 시각화 |

### 데이터 파일

| 파일 | 설명 |
|------|------|
| `01_humanness_score.csv` | position별 WT/Humanized AA, germline 일치 여부, VH3 비율 |
| `01_humanness_summary.txt` | 요약 보고서 (서열, 지표, mutation별 상세) |
| `03_md_comparison.csv` | MD 비교 통계 (실행 후 생성) |
| `03_md_comparison.txt` | MD 비교 보고서 (실행 후 생성) |

### Figures

| 파일 | 설명 |
|------|------|
| `figures/01_humanness_comparison.png/pdf` | WT vs Humanized humanness score 비교 bar plot |
| `figures/03_rmsd_comparison.png/pdf` | Backbone RMSD 시계열 비교 (실행 후 생성) |
| `figures/03_rmsf_comparison.png/pdf` | Per-residue RMSF 비교 (실행 후 생성) |
| `figures/03_rg_comparison.png/pdf` | Radius of gyration 비교 (실행 후 생성) |

---

## 의존성

```
../02_rosetta_ddg/01_imgt_pdb_mapping.csv           ← IMGT ↔ PDB mapping
../03_vh3_family_analysis/04_VH3_imgt_aa_proportion.csv  ← VH3 AA proportions
../03_vh3_family_analysis/03_VH3_unique_imgt.csv    ← IGHV3-66 IMGT residues
../01_structure_prediction/fold_humanized_fap_nb_model_0.cif ← Humanized AF3 structure
../05_gromacs_mmpbsa/*.mdp                          ← MD parameter files
```

---

## 실행

```bash
# 1. Humanness score 계산
python 01_humanness_score.py

# 2. Humanized nanobody MD 실행 (~12-24시간, GPU 필요)
bash 02_run_md_humanized.sh

# 3. xvg 파일 추출 (WT trajectory 경로 확인 필요)
bash 03_extract_xvg.sh

# 4. WT vs Humanized MD 비교
python 03_compare_md.py
```

---

## 참고 문헌

- Vincke, C., et al. (2009). General strategy to humanize a camelid single-domain antibody and identification of a universal humanized nanobody scaffold. *J. Biol. Chem.*, 284(5), 3273-3284.
- Rossotti, M.A., et al. (2022). Streamlined method for parallel identification of single domain antibodies to membrane receptors on whole cells. *Biochim. Biophys. Acta Gen. Subj.*, 1866(1), 130038.
- Prihoda, D., et al. (2022). BioPhi: A platform for antibody design, humanization, and humanness evaluation. *mAbs*, 14(1), 2020203.
