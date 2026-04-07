# VH3 Family Germline Allele Sequences

Human VH3 family (IGHV3-*) germline allele 서열 모음.
Nanobody (camelid VHH)는 human VH3 family와 구조적으로 가장 유사하여 비교 분석 목적에 활용 가능.

---

## 파일 목록

### 스크립트

| 파일 | 설명 |
|------|------|
| `01_fetch_vh3_alleles.py` | OGRDB REST API에서 VH3 allele 서열 다운로드 |
| `02_translate_to_aa.py` | DNA → 아미노산 번역 후 CSV 갱신 |
| `03_imgt_numbering.py` | unique 서열 추출 + ANARCI IMGT numbering |
| `04_aa_proportion.py` | IMGT position별 아미노산 비율 계산 |

### 데이터 파일

| 파일 | 설명 |
|------|------|
| `01_VH3_alleles_gapped.fasta` | IMGT gap 포함 DNA FASTA (93 alleles) |
| `01_VH3_alleles_ungapped.fasta` | gap 제거 DNA FASTA (93 alleles) |
| `01_VH3_alleles_summary.csv` | allele별 DNA + 아미노산 서열 (93 alleles) |
| `03_VH3_unique.fasta` | unique 아미노산 서열 FASTA (67 seqs) |
| `03_VH3_unique_imgt.csv` | unique 서열별 IMGT position 컬럼 CSV (67 seqs × 118 positions) |
| `04_VH3_imgt_aa_proportion.csv` | position별 아미노산 비율 wide format |
| `04_VH3_imgt_aa_proportion_long.csv` | position별 아미노산 비율 long format |

---

## 1. 서열 다운로드 (`01_fetch_vh3_alleles.py`)

**데이터 출처:** [OGRDB](https://ogrdb.airr-community.org/) (Open Germline Receptor Database) — AIRR Community 공식 germline DB

**API 엔드포인트:**
```
GET https://ogrdb.airr-community.org/api_v2/germline/set/9606.IGH_VDJ/latest/gapped
```

- `9606` : Homo sapiens (NCBI Taxonomy ID)
- `IGH_VDJ` : 중쇄(heavy chain) V, D, J 유전자 세트
- `latest` : 최신 버전
- `gapped` : IMGT 갭(`.`) 포함 FASTA 포맷

전체 human IGHV 서열(245개) 중 header가 `IGHV3-` 로 시작하는 것만 필터링 → **93개 allele** 추출.

**저장 형식:**
- Gapped FASTA : IMGT 갭 문자(`.`) 보존 → IMGT 넘버링 기반 정렬 분석용
- Ungapped FASTA : `.`과 `-` 제거 → 실제 coding sequence

---

## 2. 아미노산 번역 (`02_translate_to_aa.py`)

`01_VH3_alleles_summary.csv`의 `sequence` 컬럼(ungapped DNA)을 **표준 유전암호(standard genetic code)** 로 번역하여 `aa_sequence` 컬럼을 추가한다.

- 코돈(3 nt) 단위로 5'→3' 순차 번역
- 인식 불가 코돈 → `X` / 정지 코돈(TAA·TAG·TGA) → `*`
- 외부 라이브러리 미사용, 스크립트 내 `CODON_TABLE` dict로 직접 처리

예시 (`IGHV3-11*01`):
```
DNA : CAGGTGCAGCTGGTGGAGTCTGGG...
AA  : QVQLVESG...YYCAR  (98 aa)
```

---

## 3. IMGT Numbering (`03_imgt_numbering.py`)

**unique 서열 추출** → `aa_sequence` 기준 중복 제거: 93 alleles → **67 unique sequences**
- 동일 서열을 공유하는 allele는 `all_alleles` 컬럼에 `;` 구분으로 모두 보존
- 대표명(`representative`)은 첫 등장 allele

**IMGT numbering** → [ANARCI](https://github.com/oxpig/ANARCI) (Antibody Numbering And Receptor ClassIfication) 사용

```python
anarci(sequences, scheme='imgt', allow={'H'})
```

- HMMER 기반 HMM profile로 VH 도메인 인식 후 IMGT 번호 부여
- Position label: 정수 위치는 `'5'`, 삽입은 `'111A'` 형식
- gap(`-`) : 해당 서열에 그 position이 없음

출력 (`03_VH3_unique_imgt.csv`) 컬럼 구조:
```
representative | all_alleles | chain_type | aa_sequence | 1 | 2 | ... | 117
```

---

## 4. 아미노산 비율 계산 (`04_aa_proportion.py`)

`03_VH3_unique_imgt.csv`를 읽어 **IMGT position별 아미노산 출현 비율**을 계산한다.

- gap(`-`)은 proportion 계산에서 제외하고 `gap_fraction` 컬럼으로 별도 표기
- Wide format : rows = position, cols = 아미노산 (비율 0~1)
- Long format : `position`, `amino_acid`, `count`, `proportion` 4컬럼

**가장 다양한 IMGT positions (CDR2 중심):**

| IMGT position | 아미노산 종류 | top 잔기 |
|:---:|:---:|---|
| 55 | 10 | V(14), G(13), A(11) |
| 38 | 9 | A(22), Y(12), D(9) |
| 66 | 9 | Y(44), G(7), H(6) |
| 58 | 8 | S(27), W(16), Y(7) |
| 64 | 8 | S(31), T(13), N(7) |

---

## 재실행

```bash
# 1. 서열 새로 다운로드
python3 01_fetch_vh3_alleles.py

# 2. DNA → 아미노산 번역 (CSV 갱신)
python3 02_translate_to_aa.py

# 3. unique 서열 추출 + IMGT numbering
python3 03_imgt_numbering.py

# 4. position별 아미노산 비율 계산
python3 04_aa_proportion.py
```
