# 나노바디 Humanization ddG 분석 - 타임라인

## 목표
나노바디의 13개 framework 돌연변이가 안정성에 미치는 영향을 Cartesian ddG로 계산

## 최종 워크플로우

### 1. AlphaFold3 구조 확인 ✅
**시간**: 16:30
**위치**: `/home/laugh/rosetta/alphafold3_prediction_FAP_Nb/`
**방법**: 이미 예측된 AlphaFold3 구조 발견

**결과**:
- 5개 모델 발견 (model_0 ~ model_4)
- 모든 모델 동일한 confidence: **pTM = 0.88** (매우 높음!)
- Sequence length: **125 residues** (정확한 target 길이)
- Quality indicators:
  - has_clash = 0 (충돌 없음)
  - fraction_disordered = 0 (disordered region 없음)
  - 10 recycling iterations

**선택**: Model 0 사용 (`fold_anti_fap_nb_model_0.cif`)

---

### 2. IMGT 넘버링 매핑 ✅
**시간**: 16:18-16:19
**도구**: `fix_imgt_numbering.py`
**방법**: IMGT 넘버링을 sequential PDB numbering으로 변환

**Target sequence** (125 residues):
```
QVQLQESGGGSVQAGGSLRLSCAASGYTVRSSYMGWFRQVPGKQREAVAIITSGGTTYYADSVKGRF
TISRDNAKNTLYLQMNSLKPEDTAMYYCAGRTGFIGGIWFRDRDYDYWGQGTQVTVSS
```

**돌연변이 매핑 결과**:

| IMGT Position | Sequential Position | Mutation | Region | Verification |
|---------------|---------------------|----------|--------|--------------|
| 1             | 1                   | Q1E      | FR1    | ✓ MATCH      |
| 5             | 5                   | Q5V      | FR1    | ✓ MATCH      |
| 12            | 7                   | S7L      | FR1    | ✓ Found (offset -5) |
| 15            | 14                  | A14P     | FR1    | ✓ Found (offset -1) |
| 40            | 35                  | G35S     | FR2    | ✓ Found      |
| 45            | 40                  | V40A     | FR2    | ✓ Found      |
| 54            | 49                  | A49S     | FR2    | ✓ Found      |
| 55            | 50                  | I50V     | FR2    | ✓ Found      |
| 83            | 74                  | A74S     | FR3    | ✓ Found (offset -9) |
| 95            | 86                  | K86R     | FR3    | ✓ Found (offset -9) |
| 96            | 87                  | P87A     | FR3    | ✓ Found (offset -9) |
| 101           | 92                  | M92V     | FR3    | ✓ Found (offset -9) |
| 123           | 120                 | Q120L    | FR4    | ✓ MATCH      |

**결과**: 13개 돌연변이 모두 확인 완료!

**생성 파일**: `mutations_final.txt`
```
total 13
1
Q 1 E
1
Q 5 V
...
```

---

### 3. Cartesian ddG 계산 시작 ✅
**시간**: 16:50
**도구**: `cartesian_ddg.linuxgccrelease`
**입력**: AlphaFold3 model 0 (fold_anti_fap_nb_model_0.cif)

**계산 파라미터**:
```bash
cartesian_ddg.linuxgccrelease \
  -in:file:s fold_anti_fap_nb_model_0.cif \
  -ddg:mut_file mutations_final.txt \
  -ddg:iterations 5              # 5 iterations for accuracy
  -ddg:cartesian                 # Cartesian space
  -ddg:bbnbrs 1                  # Backbone flexibility ±1 residue
  -ddg:dump_pdbs false           # Don't save intermediate PDBs
  -out:path:all ddg_output_af3/
  -score:weights ref2015_cart    # Cartesian score function
  -fa_max_dis 9.0               # Repack within 9Å
  -ddg:legacy true              # Legacy mode
```

**상태**: 🔄 실행 중
- WT baseline 계산 중
- 예상 완료: 3-6시간 (총 70 calculations = 14 mutations × 5 iterations)

**출력**:
- Log: `ddg_af3_calculation.log`
- Results: `ddg_output_af3/ddg_predictions.out` (생성 중)
- Task ID: b953799

**구조 정보 (로그에서 확인)**:
- ✓ CIF 파일 정상 로드
- ✓ Disulfide bond 감지: Cys 22-95
- ✓ Score function 초기화 완료
- ✓ FastRelax 진행 중

---

## 진행 상황 요약

### ✅ 완료
1. ✅ AlphaFold3 구조 확인 및 선택 (pTM 0.88)
2. ✅ IMGT 넘버링 → Sequential 매핑 (13/13 성공)
3. ✅ 돌연변이 파일 생성 (`mutations_final.txt`)
4. ✅ Cartesian ddG 계산 시작 (AlphaFold3 구조)

### 🔄 진행 중
- 🔄 WT baseline 계산 (FastRelax)
- ⏱ 예상: WT ~20-30분, 각 돌연변이 ~15-20분
- ⏱ 총 예상 완료: 3-6시간

### ⏳ 대기
1. ⏳ ddG 계산 완료 대기
2. ⏳ 결과 분석 (`analyze_ddg_results.py`)
3. ⏳ High-risk 돌연변이 식별
4. ⏳ 실험 권장사항 도출

---

## 생성된 파일

### 입력 파일
```
nanobody.fasta                      # Target sequence
mutations_final.txt                 # 13 mutations for ddG
fold_anti_fap_nb_model_0.cif       # AlphaFold3 structure ⭐
```

### 스크립트
```
fix_imgt_numbering.py              # IMGT numbering mapper
analyze_ddg_results.py             # Results analysis
monitor_af3_ddg.sh                 # Progress monitor
calc_rmsd.py                       # Structure comparison
```

### 출력 (진행 중)
```
ddg_af3_calculation.log            # Execution log
ddg_output_af3/                    # Results directory
  └── ddg_predictions.out          # (생성 예정)
```

### 문서
```
README.md                          # Project overview
TIMELINE.md                        # This file
AF3_CALCULATION_STATUS.md          # Current status
FINAL_SUMMARY.md                   # Comprehensive summary
```

---

## 예상 결과

### High-risk 돌연변이 (ddG > +2.0 REU)
1. **A14P** (IMGT 15, position 14)
   - Proline 도입 → backbone rigidity 증가
   - Framework 1 영역
   - **권장**: Revert 강력 추천

2. **G35S** (IMGT 40, position 35)
   - Glycine flexibility 손실
   - Framework 2 영역
   - **권장**: Revert 강력 추천

3. **P87A** (IMGT 96, position 87)
   - Proline 제거
   - Framework 3 영역
   - **주의**: Cartesian ddG artifact 가능
   - **권장**: 결과 확인 후 결정

### Medium-risk 돌연변이 (ddG +1.0 ~ +2.0 REU)
- S7L, A49S, M92V: 중간 크기 변화
- **권장**: Revert 고려

### Safe 돌연변이 (ddG < +0.5 REU)
- Q1E, Q5V: Surface charges
- V40A, I50V: Similar size
- K86R: Conservative
- A74S, Q120L: Small changes
- **권장**: 유지 가능

---

## 모니터링

### 실시간 진행 확인
```bash
cd /home/laugh/rosetta/nanobody_humanization

# Progress check
./monitor_af3_ddg.sh

# Real-time log
tail -f ddg_af3_calculation.log

# Background task output
tail -f /tmp/claude-1000/-home-laugh-rosetta/tasks/b953799.output
```

### 완료 확인
```bash
# Check for completion
grep "reported success" ddg_af3_calculation.log | wc -l

# Should see 70 (14 mutations × 5 iterations)
```

---

## 다음 단계

### 1. 계산 완료 대기 (3-6시간)
- WT baseline: ~30분
- 13 mutations: ~3-5시간
- Total: ~3.5-5.5시간

### 2. 결과 분석
```bash
python3 analyze_ddg_results.py ddg_output_af3/ddg_predictions.out
```

**출력**:
- High-risk 돌연변이 목록
- ddG 값과 표준편차
- 해석 thresholds
- 실험 권장사항

### 3. 최종 권장사항
**실험 계획**:
1. **High-risk 제거 버전**: A14P, G35S 제거
2. **Medium-risk 평가**: S7L, A49S, M92V 포함/제외 버전
3. **최종 humanized**: Safe 돌연변이만 포함

---

## 중요 교훈

### 1. AlphaFold3 우선
- ✓ 항상 AlphaFold 구조 먼저 확인
- ✓ High confidence (pTM > 0.8)이면 신뢰
- ✓ Sequence length 정확성 중요

### 2. IMGT 넘버링
- ✓ IMGT와 sequential은 직접 매핑 안 됨
- ✓ 각 영역별로 offset 다름 (FR1: ~0, FR2: ~3, FR3: ~9)
- ✓ 실제 서열에서 확인 필수

### 3. Cartesian ddG
- ✓ 구조 품질이 매우 중요
- ✓ 5+ iterations 권장
- ✓ Proline 관련 artifact 주의
- ✓ Framework 돌연변이는 비교적 안정적

### 4. 계산 시간
- WT baseline: 가장 오래 걸림 (~30분)
- 각 돌연변이: ~15-20분
- Iterations 많을수록 정확하지만 시간 증가

---

**최종 업데이트**: 2026-02-09 16:51
**상태**: ✅ 완료 (Set 1)
**Task ID**: b953799

---

## Set 2: IMGT 55/96 Saturation ddG 분석

### 배경
Set 1 결과에서 두 위치가 "Strongly Destabilizing"으로 판정:
- **IMGT 55** (sequential 50, I50V): +2.52 REU
- **IMGT 96** (sequential 87, P87A): +3.44 REU

이 두 위치에 다른 아미노산을 치환하여 위치별 tolerance profile을 파악한다.

### 새 돌연변이 10개

| IMGT | Sequential | WT | 치환 |
|------|-----------|-----|------|
| 55   | 50        | I  | G, A, S, Y, R, F, L, Q |
| 96   | 87        | P  | V, D |

### 실행 (2026-02-11)

**입력 파일**: `mutations_set2.txt`
```
total 10
1
I 50 G
1
I 50 A
...
1
P 87 D
```

**명령어**:
```bash
/home/laugh/rosetta/source/bin/cartesian_ddg.linuxgccrelease \
  -in:file:s fold_anti_fap_nb_model_0.cif \
  -ddg:mut_file mutations_set2.txt \
  -ddg:iterations 5 \
  -ddg:cartesian \
  -ddg:bbnbrs 1 \
  -ddg:dump_pdbs false \
  -out:path:all ddg_output_set2/ \
  -score:weights ref2015_cart \
  -fa_max_dis 9.0 \
  -ddg:legacy true \
  > ddg_set2_calculation.log 2>&1
```

**계산 규모**: (10 mutations + 1 WT) × 5 iterations = 55 calculations
**예상 소요**: ~3시간

### 진행 상황

- ✅ mutations_set2.txt 생성
- ✅ ddg_output_set2/ 디렉토리 생성
- 🔄 cartesian_ddg 계산 중 (Task ID: b2bb26a)
- ⏳ 결과 파싱: `python3 convert_ddg_results.py mutations_set2.ddg && mv ddg_predictions.out ddg_predictions_set2.out`
- ⏳ 결과 분석: `python3 analyze_ddg_results.py ddg_predictions_set2.out`

### 완료 확인 방법
```bash
# 계산 완료 확인 (55+ lines)
grep "MUT_" ddg_set2_calculation.log | wc -l

# 결과 파일 확인
cat ddg_predictions_set2.out
```

**최종 업데이트**: 2026-02-11
**상태**: 🔄 계산 중
**Task ID**: b2bb26a
