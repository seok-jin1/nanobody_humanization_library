# AlphaFold3 기반 Cartesian ddG 계산 - 진행 상황

## 🎯 현재 상태: 실행 중 ✅

**시작 시간**: 2026-02-09 16:50
**예상 완료**: 3-6시간 후 (약 19:50 - 22:50)

## 📊 계산 설정

### 구조 정보
- **Model**: AlphaFold3 fold_anti_fap_nb_model_0.cif
- **Confidence**: pTM = 0.88 (Excellent) ⭐
- **Sequence length**: 125 residues (정확!)
- **Disulfide bonds**: Cys 22-95 (감지됨)

### 계산 파라미터
```bash
cartesian_ddg.linuxgccrelease \
  -in:file:s fold_anti_fap_nb_model_0.cif \
  -ddg:mut_file mutations_final.txt \
  -ddg:iterations 5          # 5 iterations (정확도 향상)
  -ddg:cartesian             # Cartesian space
  -ddg:bbnbrs 1              # Backbone flexibility ±1
  -out:path:all ddg_output_af3/
  -score:weights ref2015_cart
  -fa_max_dis 9.0
  -ddg:legacy true
```

### 돌연변이 목록 (13개)
```
Q1E   (position 1,   IMGT 1)    - FR1
Q5V   (position 5,   IMGT 5)    - FR1
S7L   (position 7,   IMGT 12)   - FR1
A14P  (position 14,  IMGT 15)   - FR1  ⚠️ High-risk expected
G35S  (position 35,  IMGT 40)   - FR2  ⚠️ High-risk expected
V40A  (position 40,  IMGT 45)   - FR2
A49S  (position 49,  IMGT 54)   - FR2
I50V  (position 50,  IMGT 55)   - FR2
A74S  (position 74,  IMGT 83)   - FR3
K86R  (position 86,  IMGT 95)   - FR3
P87A  (position 87,  IMGT 96)   - FR3  ⚠️ High-risk expected
M92V  (position 92,  IMGT 101)  - FR3
Q120L (position 120, IMGT 123)  - FR4
```

## 📈 진행 상황

**총 계산 수**: 14 mutations (13 + 1 WT) × 5 iterations = **70 calculations**

**현재 단계**: WT baseline 계산 중 🔄

**예상 시간**:
- WT baseline: ~20-30분
- 각 돌연변이: ~15-20분 × 13 = 3-4시간
- **총 예상**: 3.5-5시간

## 🔍 모니터링

### 실시간 확인
```bash
cd /home/laugh/rosetta/nanobody_humanization

# 진행 상황 체크
./monitor_af3_ddg.sh

# 실시간 로그
tail -f ddg_af3_calculation.log

# 백그라운드 task 출력
tail -f /tmp/claude-1000/-home-laugh-rosetta/tasks/b953799.output
```

### 완료 확인
```bash
# ddG 결과 파일 확인
ls -la ddg_output_af3/

# 결과 파일 확인
cat ddg_output_af3/ddg_predictions.out
```

## 📊 예상 결과

### High-risk 예상 돌연변이
1. **A14P** (IMGT 15)
   - Proline 도입 → backbone rigidity 증가
   - 예상 ddG: +2.0 ~ +3.5 REU

2. **G35S** (IMGT 40)
   - Glycine flexibility 손실
   - 예상 ddG: +1.5 ~ +2.5 REU

3. **P87A** (IMGT 96)
   - Proline 제거
   - 예상 ddG: +1.0 ~ +2.5 REU (artifact 가능)

### Safe 예상 돌연변이
- Q1E, Q5V: Surface charges (예상 ddG: -0.5 ~ +0.5)
- K86R: Conservative (예상 ddG: -0.2 ~ +0.5)
- V40A, I50V: Similar size (예상 ddG: 0.0 ~ +1.0)

## 🎯 완료 후 작업

### 1. 결과 분석 (자동)
```bash
python3 analyze_ddg_results.py ddg_output_af3/ddg_predictions.out > af3_analysis_report.txt
```

### 2. 결과 해석
**기준값:**
- **ddG > +2.0 REU**: ❌ 강력히 revert 추천
- **ddG +1.0 ~ +2.0 REU**: ⚠️ Revert 고려
- **ddG +0.5 ~ +1.0 REU**: ⚡ 모니터링
- **ddG < +0.5 REU**: ✅ 안전

### 3. 최종 권장사항 생성
실험 우선순위:
1. High-risk 돌연변이 제거 버전 테스트
2. Medium-risk는 유지하고 발현량 확인
3. Safe 돌연변이만 포함한 최종 humanized 버전

## 📁 파일 구조

```
nanobody_humanization/
├── AF3_CALCULATION_STATUS.md       # 이 파일
├── fold_anti_fap_nb_model_0.cif    # AlphaFold3 구조
├── mutations_final.txt             # 돌연변이 리스트
├── ddg_af3_calculation.log         # 실행 로그
├── monitor_af3_ddg.sh              # 모니터링 스크립트
├── analyze_ddg_results.py          # 분석 스크립트
├── ddg_output_af3/                 # 결과 디렉토리 (생성 중)
│   └── ddg_predictions.out         # (계산 완료 후)
└── archive/
    └── template_based/             # 1ZVH 기반 계산 (참고용)
```

## 🔬 왜 AlphaFold3인가?

### Template (1ZVH) vs AlphaFold3

| 특성 | Template | AlphaFold3 | 승자 |
|------|----------|------------|------|
| Sequence length | 127 ❌ | 125 ✓ | AF3 |
| Confidence | N/A | pTM 0.88 | AF3 |
| RMSD to AF3 | 18.13 Å | 0.0 Å | AF3 |
| CDR quality | Moderate | Excellent | AF3 |
| Framework quality | Good | Excellent | AF3 |
| Bias | Template bias | De novo | AF3 |
| **권장도** | 참고용 | **메인** ⭐ | **AF3** |

### 핵심 이유
1. ✓ **정확한 길이**: 125 residues (template은 127)
2. ✓ **높은 신뢰도**: pTM 0.88
3. ✓ **구조 차이**: Template과 18 Å RMSD → 다른 예측 가능
4. ✓ **No bias**: De novo prediction

## 📖 참고 문서

- `FINAL_SUMMARY.md` - 전체 프로젝트 종합
- `TIMELINE.md` - 시도한 방법들
- `README.md` - 프로젝트 개요
- `calc_rmsd.py` - 구조 비교 결과

## ⚠️ 중요 노트

1. **계산 시간**: 백그라운드로 실행 중, 중단하지 마세요!
2. **결과 신뢰도**: AlphaFold3 구조 사용으로 template보다 정확
3. **Proline artifact**: P87A는 여전히 artifact 가능성 있음
4. **통계적 신뢰도**: 5 iterations로 평균과 표준편차 계산

---

**마지막 업데이트**: 2026-02-09 16:51
**상태**: WT baseline 계산 중 🔄
**예상 완료**: ~3-6시간
