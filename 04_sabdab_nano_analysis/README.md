# SAbDab Nanobody Structure Analysis

약 1,900여 개의 nanobody 구조(SAbDab-nano)를 대상으로 IMGT 넘버링 기반 잔기 분포, 변이도(Shannon entropy), 아미노산 물성 비율, CDR 길이 분포를 분석한다.

---

## 데이터 준비

### Input 데이터

| 파일 | 설명 |
|------|------|
| `sabdab_nano_summary_all.tsv` | SAbDab nanobody summary (PDB ID, Hchain 매핑 등) |
| `imgt/*.pdb` | IMGT 리넘버링된 PDB 구조 파일 (~1,900개) |

`imgt/` 디렉토리는 용량 문제로 repository에 포함하지 않는다.
원본 위치: `/home/laugh/nanobody_residue/imgt` (SAbDab에서 다운로드 후 IMGT 리넘버링)

---

## 파일 목록

### 스크립트

| 파일 | 설명 |
|------|------|
| `01_analyze_and_plot_imgt.py` | IMGT PDB 파싱 → 잔기 빈도 매트릭스 생성 + 히트맵 시각화 |
| `02_create_original_summary.py` | 원본 IMGT position별 빈도 20% 이상 잔기 요약 |
| `03_create_numeric_summary.py` | 숫자 position 그룹화 후 빈도 20% 이상 잔기 요약 |
| `04_visualize_numeric_only.py` | Shannon entropy + 아미노산 물성 비율 시각화 |
| `05_analyze_cdr3_length.py` | CDR1/CDR2/CDR3 길이 분포 계산 + 시각화 |

### 데이터 파일

| 파일 | 설명 |
|------|------|
| `01_imgt_residue_matrix.csv` | 위치별 아미노산 빈도 비율 매트릭스 |
| `02_imgt_original_summary_top20.csv` | 원본 position 기준 주요 잔기 요약 (>=20%) |
| `03_imgt_numeric_summary_top20.csv` | 숫자 position 기준 주요 잔기 요약 (>=20%) |
| `04_plot_data.csv` | entropy + 물성 비율 plot data |
| `05_cdr_length_distribution.csv` | 각 nanobody의 CDR1/2/3 길이 및 서열 |

### Figures

| 파일 | 설명 |
|------|------|
| `figures/01_imgt_residue_heatmap.png/pdf` | V-domain 잔기 분포 히트맵 (IMGT 1-128) |
| `figures/04_imgt_numeric_only_analysis.png/pdf` | Shannon entropy + 물성 비율 2-panel plot |
| `figures/05_cdr_length_distribution.png/pdf` | CDR1/2/3 길이 분포 히스토그램 (3-panel) |

### 참고 문서

| 파일 | 설명 |
|------|------|
| `IMGT_numbering_ranges.md` | IMGT 기준 FR1-4, CDR1-3 영역 번호 범위 정의 |

---

## 실행 순서

```bash
# 1. IMGT PDB 파싱 → 잔기 빈도 매트릭스 + 히트맵
#    (imgt/ 디렉토리 필요)
python 01_analyze_and_plot_imgt.py

# 2. 원본 position 기준 주요 잔기 요약
python 02_create_original_summary.py

# 3. 숫자 position 기준 주요 잔기 요약
python 03_create_numeric_summary.py

# 4. Shannon entropy + 물성 비율 시각화
python 04_visualize_numeric_only.py

# 5. CDR 길이 분포 계산 + 시각화
#    (imgt/ 디렉토리 필요)
python 05_analyze_cdr3_length.py
```

01, 05는 `imgt/` PDB 디렉토리가 필요하며, 02-04는 `01_imgt_residue_matrix.csv`만 있으면 실행 가능하다.

---

## 주요 분석 결과

### Shannon Entropy
- CDR 영역 평균 entropy >> FR 영역 평균 entropy
- CDR3가 가장 높은 서열 다양성을 보임

### CDR 길이 분포
- CDR1: 대부분 8 residues (IMGT 27-38)
- CDR2: 대부분 8 residues (IMGT 56-65, insertion 포함 시 변동)
- CDR3: 14-21 residues로 길이 다양성이 가장 높음 (IMGT 105-117, insertion code 포함)
