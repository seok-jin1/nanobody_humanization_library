# Rosetta WT vs Humanized Stability Comparison

인간화 전후 nanobody의 구조적 안정성을 Rosetta로 비교 분석한다.

---

## 분석 내용

### 1. Multi-Mutant Cartesian ddG

WT 구조에 8개 humanization mutation을 **동시 적용**하여 전체 인간화의 ddG를 계산.
개별 mutation ddG 합산과 비교하면 epistatic effect(잔기간 상호작용) 확인 가능.

### 2. FastRelax Total Score 비교

WT와 Humanized AF3 구조를 각각 Rosetta FastRelax (ref2015_cart)로 에너지 최소화 후
total_score를 비교. 낮을수록 안정한 구조.

---

## 파일 목록

### 스크립트

| 파일 | 설명 |
|------|------|
| `01_run_multi_mutant_ddg.sh` | 8-mutation 동시 적용 Cartesian ddG 실행 |
| `02_run_relax_compare.sh` | WT/Humanized FastRelax 실행 (각 5 구조) |
| `03_analyze_relax.py` | 결과 파싱, 비교 보고서, box plot 시각화 |

### Input

| 파일 | 설명 |
|------|------|
| `01_multi_mutant_ddg.txt` | 8개 mutation Rosetta format |

### Output (실행 후 생성)

| 파일 | 설명 |
|------|------|
| `01_multi_mutant.ddg` | multi-mutant ddG raw output |
| `02_relax_wt.sc` | WT relaxed score (5 structures) |
| `02_relax_humanized.sc` | Humanized relaxed score (5 structures) |
| `03_stability_comparison.csv` | 비교 데이터 |
| `03_stability_comparison.txt` | 텍스트 보고서 |
| `figures/03_stability_comparison.png/pdf` | total score box plot |

---

## 실행

```bash
# 1. Multi-mutant ddG (~30분)
bash 01_run_multi_mutant_ddg.sh

# 2. FastRelax WT + Humanized (~1시간)
bash 02_run_relax_compare.sh

# 3. 결과 분석 + 시각화
python 03_analyze_relax.py
```
