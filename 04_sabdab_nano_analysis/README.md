# Nanobody Residue Analysis Project

이 프로젝트는 약 1,900여 개의 Nanobody 구조 데이터를 바탕으로 IMGT 넘버링 기준의 잔기 분포, 변이도(Entropy), 및 농축도(Enrichment)를 분석한 프로젝트입니다.

## 진행 순서 (Workflows)

### 1. 데이터 준비 및 구조화
- `all_nano_structures.zip` 압축 해제 및 데이터 확인.
- 원본 PDB 파일을 `pdbs/` 폴더로 정리.
- `sabdab_nano_summary_all.tsv`를 통해 각 PDB ID별 Nanobody 체인(Hchain) 매핑 정보 구축.

### 2. 잔기 비율 매트릭스 생성 (`analyze_and_plot_imgt.py`)
- IMGT 리넘버링된 PDB 파일을 파싱하여 위치별 아미노산 빈도 계산.
- Nanobody V-domain의 표준 범위인 **IMGT 1-128**번에 집중하여 분석.
- 결과물: `imgt_residue_matrix.csv`, `imgt_residue_heatmap.png` (전체 히트맵).

### 3. 통계적 요약 데이터 추출
- **20% 이상 빈도 잔기 추출:** 각 위치에서 지배적인 잔기들을 요약.
  - `imgt_numeric_summary_top20.csv`: 숫자형 포지션 기준 통합 요약.
  - `imgt_original_summary_top20.csv`: Insertion code(A, B, C...)를 포함한 원본 포지션 기준 요약.

### 4. 고급 시각화 분석
- **변이도 및 성질 분석 (`visualize_advanced.py`):**
  - **Shannon Entropy:** 각 위치의 서열 다양성 측정 (CDR 영역에서 높게 나타남).
  - **Physicochemical Properties:** 아미노산을 화학적 성질(疎水性, 極性 등)별로 묶어 분포 시각화.
  - 결과물: `imgt_advanced_analysis.png`.
- **숫자형 통합 시각화 (`visualize_numeric_only.py`):**
  - 삽입 코드를 제거하고 숫자 단위로 데이터를 합산하여 전체적인 경향성을 파악하기 쉽게 개선.
  - 결과물: `imgt_numeric_only_analysis.png`.

### 5. 농축도(Enrichment) 분석 (`visualize_enrichment.py`)
- **Log-Odds 분석:** Nanobody 전체 평균 아미노산 빈도(Background) 대비 특정 위치에서 유의미하게 많이 나타나는 잔기(Signature)를 발굴.
- Log2 Fold Change를 계산하여 시각화.
- 결과물: `imgt_enrichment_analysis.png`.

### 6. CDR 길이 분포 분석 (`analyze_cdr3_length.py`)
- IMGT 기준 CDR1 (27–38), CDR2 (56–65), CDR3 (105–117, 삽입 코드 포함)의 잔기 수를 각 nanobody 구조별로 계산.
- 결과물:
  - `cdr_length_distribution.csv`: 각 구조의 CDR1/2/3 길이 및 서열 (`pdb_id`, `chain`, `cdr1_length`, `cdr1_sequence`, `cdr2_length`, `cdr2_sequence`, `cdr3_length`, `cdr3_sequence` 컬럼).
  - `cdr_length_distribution.png`: CDR1/2/3 길이 빈도 분포 히스토그램 (3-panel).

### 7. 기준 문서 작성
- `IMGT_numbering_ranges.md`: 분석에 사용된 IMGT 기준 FR1-4 및 CDR1-3 영역의 번호 범위 정의.

### 8. IMGT 위치별 Shannon Entropy 계산 및 매트릭스 통합 (스크립트 없음, 인라인 처리)
- **Input:** `imgt_residue_matrix.csv`
- 각 IMGT 위치(열)에 대해 Shannon Entropy `H = -Σ p·log₂(p)` 계산 (빈도 0인 아미노산 제외).
- 기존 잔기 비율 매트릭스에 `Entropy` 행을 추가하여 통합.
- **Output:** `imgt_residue_matrix_with_entropy.csv` (24행 × 576열: 기존 23 아미노산 + Entropy 행)

### 9. IMGT 위치별 Entropy 히스토그램 시각화 (스크립트 없음, 인라인 처리)
- **Input:** `imgt_residue_matrix_with_entropy.csv`
- IMGT 1–128번 위치를 대상으로 Entropy 분포를 3-panel 히스토그램으로 시각화.
  - Panel 1: 전체 위치 (n=206, insertion code 포함)
  - Panel 2: CDR 영역만 (CDR1: 27–38, CDR2: 56–65, CDR3: 105–117, n=69)
  - Panel 3: Framework(FR) 영역만 (FR1–FR4, n=137)
- CDR 평균 Entropy (1.54 bits) >> FR 평균 Entropy (0.54 bits) → CDR 영역의 높은 서열 다양성 확인.
- **Output:** `imgt_entropy_histogram.png`

## 주요 분석 결과물
- `imgt_residue_matrix.csv`: 전체 위치별 아미노산 비율 데이터.
- `imgt_residue_matrix_with_entropy.csv`: 아미노산 비율 + IMGT 위치별 Shannon Entropy 통합 데이터.
- `imgt_numeric_summary_top20.csv`: 주요 잔기 요약본.
- `imgt_advanced_analysis.png`: 엔트로피 및 화학적 성질 분포 그래프.
- `imgt_enrichment_analysis.png`: 위치별 특이 잔기 농축 히트맵.
- `imgt_entropy_histogram.png`: IMGT 1–128 위치의 Entropy 분포 히스토그램 (전체/CDR/FR 3-panel).
- `cdr_length_distribution.csv`: 각 nanobody 구조의 CDR1/2/3 길이 및 서열.
- `cdr_length_distribution.png`: CDR1/2/3 길이 분포 히스토그램 (3-panel).
