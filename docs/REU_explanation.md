# Rosetta Energy Units (REU) 설명

## REU란?

**REU (Rosetta Energy Units)** = Rosetta의 자체 에너지 단위
- 무차원 단위 (dimensionless)
- kcal/mol이 **아님**
- Rosetta score function의 가중치 조합으로 계산됨

## kcal/mol로 변환?

### 대략적인 변환
문헌에 따르면:
- **1 REU ≈ 0.5-1.0 kcal/mol** (매우 대략적)
- 하지만 이 변환은 **신뢰도가 낮음**

### 왜 직접 변환이 어려운가?
1. REU는 여러 에너지 항의 가중합
   - van der Waals (fa_atr, fa_rep)
   - Solvation (fa_sol)
   - Hydrogen bonds (hbond_*)
   - Electrostatics (fa_elec)
   - Torsional angles (rama, omega)
   - etc.

2. 각 항마다 임의의 가중치가 적용됨
3. 실험적 ddG와의 상관관계는 있지만 1:1 대응은 아님

## 실험적 ddG와의 상관관계

### Cartesian ddG 성능
- **상관계수 (R):** ~0.6-0.7
- **RMSE:** ~1.5-2.0 kcal/mol (실험값 대비)

이는:
- Rosetta ddG가 **경향성은 잘 예측**
- **절대값은 정확하지 않음**

### 해석 기준

| Rosetta ddG (REU) | 실험적 의미 | 실제 대략 kcal/mol |
|-------------------|-------------|-------------------|
| **> +2.0** | 강한 불안정화 | > +1-2 kcal/mol |
| **+1.0 ~ +2.0** | 중간 불안정화 | +0.5-1.5 kcal/mol |
| **+0.5 ~ +1.0** | 약한 불안정화 | +0.2-0.8 kcal/mol |
| **-0.5 ~ +0.5** | 중립적 | ±0.5 kcal/mol |
| **< -0.5** | 안정화 | < -0.3 kcal/mol |

## 우리 결과 해석

### High-risk 돌연변이
```
A14P:  +6.04 REU  →  실험적으로 ~3-6 kcal/mol 불안정화 예상
A49S:  +4.16 REU  →  실험적으로 ~2-4 kcal/mol 불안정화 예상
P87A:  +3.44 REU  →  실험적으로 ~1.5-3 kcal/mol 불안정화 예상
```

→ **발현량 감소 또는 aggregation 위험**

### Safe 돌연변이
```
S11L:  -3.54 REU  →  실험적으로 ~1-3 kcal/mol 안정화 예상
G35S:  -1.38 REU  →  실험적으로 ~0.5-1.5 kcal/mol 안정화 예상
```

→ **안정성 향상 또는 중립적**

## 중요한 점

### ✓ REU로 할 수 있는 것
1. **상대적 비교** - 어느 돌연변이가 더 나쁜지/좋은지
2. **경향성 예측** - 불안정화 vs 안정화
3. **우선순위 결정** - 어느 것을 먼저 제거할지

### ✗ REU로 할 수 없는 것
1. **정확한 실험값 예측** - 절대값은 부정확
2. **직접 열역학 계산** - Tm, Kd 등은 예측 불가
3. **kcal/mol로 정확한 변환** - 1:1 대응 없음

## 실제 실험 검증 필요

Rosetta ddG는:
- **스크리닝 도구**로 최적
- **정확한 정량적 예측은 불가**
- **실험적 검증 필수**

### 실험 계획
1. High-risk 돌연변이 제거 버전 제작
2. 발현량, 안정성 (Tm, aggregation) 측정
3. Rosetta 예측과 비교

## 문헌 참고

- Kellogg et al. (2011): "Role of conformational sampling in computing mutation-induced changes in protein structure and stability"
- Park et al. (2016): "Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules"

**결론:** REU는 kcal/mol이 아니며, 직접 변환은 신뢰도가 낮습니다. 
상대적 비교와 경향성 파악에 사용하고, 정량적 예측은 실험으로 확인해야 합니다.
