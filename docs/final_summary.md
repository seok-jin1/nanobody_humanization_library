# Nanobody Humanization - ddG 분석 최종 결과

## 계산 정보
- ✅ 총 13개 humanization 돌연변이 분석 완료
- ✅ S11L (MUT_11LEU) 정확히 계산됨
- ✅ 각 돌연변이 5 iterations (총 70 calculations)

## 전체 요약

| 카테고리 | 개수 | 비율 |
|---------|------|------|
| **강한 불안정화** (ddG > 2.0) | 6 | 46.2% |
| **중립** (-0.5 < ddG < 0.5) | 2 | 15.4% |
| **안정화** (ddG < -0.5) | 5 | 38.5% |

---

## 🚨 HIGH RISK - 제거 강력 권장 (6개)

| 순위 | Mutation | IMGT | Region | ddG (REU) | Std | 평가 |
|------|----------|------|--------|-----------|-----|------|
| 1 | **A14P** | 15 | FR1 | **+6.04** | 0.37 | ⚠️ Proline 도입 - backbone 제약 |
| 2 | **A49S** | 54 | FR2 | **+4.16** | 1.25 | ⚠️ 강한 불안정화 |
| 3 | **P87A** | 96 | FR3 | **+3.44** | 0.91 | ⚠️ Proline 제거 (artifact 가능) |
| 4 | **V40A** | 45 | FR2 | **+2.76** | 0.17 | ⚠️ 중간 불안정화 |
| 5 | **I50V** | 55 | FR2 | **+2.52** | 0.20 | ⚠️ 중간 불안정화 |
| 6 | **Q1E** | 1 | FR1 | **+2.09** | 0.31 | ⚠️ N-terminus 위치 |

---

## ✓ NEUTRAL - 영향 적음 (2개)

| Mutation | IMGT | Region | ddG (REU) | Std | 평가 |
|----------|------|--------|-----------|-----|------|
| K86R | 95 | FR3 | +0.23 | 1.80 | Conservative 치환 |
| Q5V | 5 | FR1 | -0.09 | 0.31 | 중립적 |

---

## ✅ SAFE - 안전/안정화 (5개)

| Mutation | IMGT | Region | ddG (REU) | Std | 평가 |
|----------|------|--------|-----------|-----|------|
| Q120L | 123 | FR4 | -0.59 | 0.91 | 약간 안정화 |
| **G35S** | 40 | FR2 | **-1.38** | 0.36 | 안정화 ✓ |
| A74S | 83 | FR3 | -2.13 | 0.23 | 안정화 ✓ |
| M92V | 101 | FR3 | -2.27 | 0.83 | 안정화 ✓ |
| **S11L** | 12 | FR1 | **-3.54** | 2.95 | 강한 안정화 ✓ |

---

## 🔍 주요 발견

### 1. **A14P (IMGT 15) - 최악의 돌연변이**
- ddG = +6.04 REU (가장 높음)
- Proline 도입으로 backbone flexibility 크게 제약
- **반드시 제거 권장**

### 2. **S11L (IMGT 12) - 최고의 돌연변이**
- ddG = -3.54 REU (가장 낮음)
- 강한 안정화 효과
- **반드시 유지 권장**
- 참고: 이전 잘못된 S7L에서 S11L로 수정됨 ✓

### 3. **P87A (IMGT 96) - 주의 필요**
- ddG = +3.44 REU
- Proline 제거는 Cartesian ddG에서 artifact 가능성
- 실험적 검증 필요

### 4. **G35S (IMGT 40) - 예상과 다른 결과**
- ddG = -1.38 REU (안정화!)
- 예상: Glycine flexibility 손실로 불안정화
- 실제: 안정화 효과
- Serine이 더 나은 local structure 형성 가능

### 5. **FR2 region - 문제 집중**
- FR2의 4개 돌연변이 중 3개가 high risk
- V40A (+2.76), A49S (+4.16), I50V (+2.52)
- FR2 구조가 특히 민감한 것으로 판단

---

## 📋 실험 권장사항

### Phase 1: High-risk 제거 버전
제거할 돌연변이:
- ❌ A14P (IMGT 15)
- ❌ A49S (IMGT 54)
- ❌ P87A (IMGT 96)
- ❌ V40A (IMGT 45)
- ❌ I50V (IMGT 55)
- ❌ Q1E (IMGT 1)

유지할 돌연변이 (7개):
- ✅ Q5V (IMGT 5)
- ✅ **S11L (IMGT 12)** ⭐
- ✅ G35S (IMGT 40)
- ✅ A74S (IMGT 83)
- ✅ K86R (IMGT 95)
- ✅ M92V (IMGT 101)
- ✅ Q120L (IMGT 123)

### Phase 2: P87A 검증
- P87A는 artifact 가능성 있으므로 별도 실험 필요
- P87A 포함/제외 버전 비교

### Phase 3: 최종 최적화
- Phase 1, 2 결과 기반으로 최종 humanized 버전 확정

---

## ⚠️ 중요 주의사항

1. **표준편차가 큰 돌연변이** (std > 1.0):
   - MUT_49SER (1.25)
   - MUT_86ARG (1.80)
   - **MUT_11LEU (2.95)** - 더 많은 iteration 필요 가능
   
2. **Proline 관련**:
   - A14P: Proline 도입 → 진짜 불안정화
   - P87A: Proline 제거 → artifact 가능성

3. **구조적 context**:
   - Surface mutations은 더 관대
   - Buried mutations은 더 민감
   - 각 돌연변이의 구조적 위치 확인 필요

---

## 📊 통계 요약

```
WT baseline: 0.00 ± 0.00 REU (정상)

High Risk (ddG > 2.0):  6 mutations (46%)
Medium Risk (1-2):      0 mutations (0%)
Low Risk (0.5-1):       0 mutations (0%)
Neutral (-0.5 to 0.5):  2 mutations (15%)
Stabilizing (< -0.5):   5 mutations (38%)
```

**결론:** 13개 humanization 돌연변이 중 7개는 안전하지만, 6개는 높은 불안정화를 야기합니다. 특히 A14P, A49S, P87A는 발현량이나 안정성에 심각한 영향을 줄 수 있습니다.

---

## 🔬 Set 2: IMGT 55/96 포지션 Saturation 분석 (진행 중)

**목적**: Set 1에서 "Strongly Destabilizing"으로 판정된 I50V (IMGT 55, +2.52 REU)와 P87A (IMGT 96, +3.44 REU)의 위치별 tolerance profile 파악

### 새 돌연변이 목록 (10개)

| IMGT | Sequential | WT | 치환 | 기존 결과 |
|------|-----------|-----|------|-----------|
| 55   | 50        | I  | G, A, S, Y, R, F, L, Q | I50V = +2.52 REU |
| 96   | 87        | P  | V, D | P87A = +3.44 REU |

### 계산 파라미터 (Set 1과 동일)
- Iterations: 5 (총 55 calculations = 11 × 5)
- 구조: fold_anti_fap_nb_model_0.cif
- Score function: ref2015_cart

### 상태
- 🔄 **ddg_set2_calculation.log** 생성 중
- 🔄 mutations_set2.txt → ddg_output_set2/
- ⏳ 완료 후: `convert_ddg_results.py mutations_set2.ddg` → `ddg_predictions_set2.out`
- ⏳ 완료 후: `analyze_ddg_results.py ddg_predictions_set2.out`

### 예상 해석 방향
- **IMGT 55 (I50)**: 소수성 Ile 위치 - Gly/Ala/Ser(소형극성)은 불안정, Phe/Leu(유사소수성)은 양호 예상
- **IMGT 96 (P87)**: Pro 위치 - Val/Asp 모두 Pro 제거이므로 높은 ddG 예상 (P87A처럼 artifact 가능성 검토 필요)

---

**분석 일자:** 2026-02-10 (Set 1), 2026-02-11 (Set 2 진행 중)
**계산 방법:** Rosetta Cartesian ddG (legacy mode)
**구조:** AlphaFold3 model 0 (pTM 0.88)
**Iterations:** 5 per mutation
