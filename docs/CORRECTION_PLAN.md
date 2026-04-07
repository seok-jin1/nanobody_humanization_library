# S7L → S11L 수정 계획

## 문제 확인

### 현재 mutations_final.txt (잘못됨)
```
IMGT 12 → Sequential 7: S7L  ❌ 잘못됨!
```

### 정확한 매핑
```
IMGT 12 → Sequential 11: S11L  ✓ 올바름!
```

### 서열 확인
```
Position:  1234567890123456789012...
Sequence:  QVQLQESGGGSVQAGGSLRLSC...
                     ^
                     11번 위치 = S (Serine)
```

---

## 수정 계획

### Phase 1: 현재 계산 완료 대기 ✅
- 현재 실행 중인 ddG 계산 완료까지 대기
- 예상 완료: ~19:20
- Task ID: b953799

### Phase 2: S11L 단독 계산 준비

**1. 새 mutations 파일 생성**
```bash
cat > mutation_S11L_only.txt << 'EOF'
total 1
1
S 11 L
EOF
```

**2. S11L 계산 실행**
```bash
cartesian_ddg.linuxgccrelease \
  -in:file:s fold_anti_fap_nb_model_0.cif \
  -ddg:mut_file mutation_S11L_only.txt \
  -ddg:iterations 5 \
  -ddg:cartesian \
  -ddg:bbnbrs 1 \
  -ddg:dump_pdbs false \
  -out:path:all ddg_output_S11L/ \
  -score:weights ref2015_cart \
  -fa_max_dis 9.0 \
  -ddg:legacy true
```

**예상 시간:** ~10-15분

### Phase 3: 결과 교체

**1. 원본 결과 백업**
```bash
cp ddg_output_af3/ddg_predictions.out ddg_output_af3/ddg_predictions.out.backup
```

**2. S7L 결과를 S11L 결과로 교체**

**방법 A: 수동 교체 (권장)**
```bash
# S7L 결과 제거
grep -v "S7L\|7 LEU" ddg_output_af3/ddg_predictions.out > temp.out

# S11L 결과 추가
cat ddg_output_S11L/ddg_predictions.out | grep "S11L\|11 LEU" >> temp.out

# 정렬 및 교체
sort temp.out > ddg_output_af3/ddg_predictions_corrected.out
```

**방법 B: Python 스크립트 (자동화)**
```python
# merge_ddg_results.py 생성 필요
```

**3. 최종 파일 확인**
```bash
# S11L이 포함되었는지 확인
grep "11 LEU\|S11L" ddg_output_af3/ddg_predictions_corrected.out

# S7L이 제거되었는지 확인
grep "7 LEU\|S7L" ddg_output_af3/ddg_predictions_corrected.out
# (아무것도 나오지 않아야 함)
```

### Phase 4: 결과 재분석

**최종 분석 실행**
```bash
python3 analyze_ddg_results.py ddg_output_af3/ddg_predictions_corrected.out > final_analysis.txt
```

---

## 수정된 돌연변이 리스트 (최종)

| IMGT | Seq | Mutation | Region | Notes |
|------|-----|----------|--------|-------|
| 1    | 1   | Q1E      | FR1    | ✓ Correct |
| 5    | 5   | Q5V      | FR1    | ✓ Correct |
| 12   | **11** | **S11L** | FR1    | ✓ **CORRECTED** |
| 15   | 14  | A14P     | FR1    | ✓ Correct |
| 40   | 35  | G35S     | FR2    | ✓ Correct |
| 45   | 40  | V40A     | FR2    | ✓ Correct |
| 54   | 49  | A49S     | FR2    | ✓ Correct |
| 55   | 50  | I50V     | FR2    | ✓ Correct |
| 83   | 74  | A74S     | FR3    | ✓ Correct |
| 95   | 86  | K86R     | FR3    | ✓ Correct |
| 96   | 87  | P87A     | FR3    | ✓ Correct |
| 101  | 92  | M92V     | FR3    | ✓ Correct |
| 123  | 120 | Q120L    | FR4    | ✓ Correct |

---

## 실행 체크리스트

### ✅ 현재 계산 완료 후

- [ ] 현재 ddG 계산 완료 확인
- [ ] 결과 파일 백업
- [ ] S7L 결과 확인 (참고용으로 보관)

### 🔄 S11L 계산

- [ ] mutation_S11L_only.txt 생성
- [ ] S11L ddG 계산 실행
- [ ] 결과 확인

### 🔀 결과 교체

- [ ] 원본 백업
- [ ] S7L 제거
- [ ] S11L 추가
- [ ] 최종 파일 검증

### 📊 최종 분석

- [ ] 수정된 결과로 재분석
- [ ] 최종 리포트 생성
- [ ] TIMELINE.md 업데이트

---

## 자동화 스크립트 (옵션)

### merge_ddg_results.py
```python
#!/usr/bin/env python3
"""
Merge S11L results into the main ddG output, replacing S7L.
"""

import sys

def merge_results(original_file, s11l_file, output_file):
    """Merge S11L results, removing S7L."""

    # Read original results, excluding S7L
    with open(original_file, 'r') as f:
        original_lines = [line for line in f
                         if 'S7L' not in line and '7 LEU' not in line]

    # Read S11L results
    with open(s11l_file, 'r') as f:
        s11l_lines = [line for line in f
                     if 'S11L' in line or '11 LEU' in line]

    # Merge
    all_lines = original_lines + s11l_lines
    all_lines.sort()

    # Write output
    with open(output_file, 'w') as f:
        f.writelines(all_lines)

    print(f"✓ Merged results written to: {output_file}")
    print(f"  - Removed S7L entries")
    print(f"  - Added S11L entries: {len(s11l_lines)}")

if __name__ == "__main__":
    merge_results(
        "ddg_output_af3/ddg_predictions.out",
        "ddg_output_S11L/ddg_predictions.out",
        "ddg_output_af3/ddg_predictions_corrected.out"
    )
```

---

## 예상 타임라인

```
현재 (17:23) → 계산 완료 (~19:20) → S11L 계산 (~19:35) → 결과 교체 (19:40) → 완료!
          └─ 2시간 대기 ─┘    └─ 15분 ─┘         └─ 5분 ─┘

총 예상: 2시간 20분
```

---

## 중요 노트

### 왜 S7L이 선택되었나?
```python
# fix_imgt_numbering.py에서:
# IMGT 12 → Sequential 12로 먼저 체크
# Position 12 = 'V' (mismatch!)
# → 주변에서 'S' 찾기
# → Position 7에서 'S' 발견
# → 잘못된 매핑!

# 올바른 방법:
# IMGT 12는 FR1 끝부분
# FR1은 일반적으로 offset이 작음
# Position 11이 정확함 (IMGT 12 - 1 = Sequential 11)
```

### 영향 분석
- S7L과 S11L은 다른 위치
- Position 7: GGGSV... 의 S (FR1 초반)
- Position 11: ...GGGSVQA... 의 S (FR1 중반)
- 구조적 context가 다를 수 있음
- ddG 값도 다를 것으로 예상

---

**작성**: 2026-02-09 17:25
**상태**: 계획 수립 완료, 실행 대기 중
**담당**: 현재 계산 완료 후 진행
