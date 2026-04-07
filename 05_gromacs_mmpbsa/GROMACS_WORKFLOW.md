# GROMACS → gmx_MMPBSA Workflow
## Nanobody Mutation Stability Analysis

---

## Overview

Single-trajectory MMPBSA approach:
1. WT 단백질로 100 ns MD 1회 실행
2. 동일 trajectory에서 13개 mutation in silico 적용
3. ΔΔG = ΔG(mutant) - ΔG(WT) in **kcal/mol**

**실행 환경:**
- GROMACS 2024.1 (`/usr/local/gromacs/bin/gmx`) — CUDA GPU 지원
- gmx_MMPBSA 1.6.4 (`conda activate docking-md`)
- GPU: RTX 4060 Ti

---

## Step 1: 구조 변환 (CIF → PDB)

**입력:** `fold_anti_fap_nb_model_0.cif` (AlphaFold3, pTM=0.88)

```python
# BioPython MMCIFParser 사용
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('nb', 'fold_anti_fap_nb_model_0.cif')
io = PDBIO()
io.set_structure(structure)
io.save('gromacs_mmpbsa/nanobody_af3.pdb')
```

**출력:** `nanobody_af3.pdb`
- Chain A, 125 residues, 965 atoms

---

## Step 2: 토폴로지 생성 (`pdb2gmx`)

**입력:** `nanobody_af3.pdb`

```bash
GMX=/usr/local/gromacs/bin/gmx
cd gromacs_mmpbsa/

$GMX pdb2gmx \
  -f nanobody_af3.pdb \
  -o nanobody_processed.gro \
  -p topol.top \
  -i posre.itp \
  -ignh \
  -ff amber99sb-ildn \
  -water tip3p
```

| 파라미터 | 값 | 설명 |
|---------|-----|------|
| `-ff` | `amber99sb-ildn` | 단백질 MD 표준 force field |
| `-water` | `tip3p` | 물 모델 |
| `-ignh` | — | 기존 H 무시, 자동 추가 |

**출력:**
- `nanobody_processed.gro` — 좌표 파일 (H 추가, 1888 atoms)
- `topol.top` — 토폴로지
- `posre.itp` — position restraints

**주요 결과:**
- Cys22-Cys95 disulfide bond 자동 감지 및 링크 ✅
- Total charge: **+4e**

---

## Step 3: 시뮬레이션 박스 생성 (`editconf`)

**입력:** `nanobody_processed.gro`

```bash
$GMX editconf \
  -f nanobody_processed.gro \
  -o nanobody_box.gro \
  -c -d 1.2 -bt dodecahedron
```

| 파라미터 | 값 | 설명 |
|---------|-----|------|
| `-d` | `1.2` | 단백질 표면에서 최소 1.2 nm |
| `-bt` | `dodecahedron` | 최소 이미지 박스 (물 분자 절약) |

**출력:** `nanobody_box.gro`
- Box volume: **245.60 nm³**

---

## Step 4: 용매화 (`solvate`)

**입력:** `nanobody_box.gro`

```bash
$GMX solvate \
  -cp nanobody_box.gro \
  -cs spc216.gro \
  -o nanobody_solv.gro \
  -p topol.top
```

**출력:** `nanobody_solv.gro`
- **7,442개** TIP3P 물 분자 추가
- 밀도: 1005.81 g/L

---

## Step 5: 이온 추가 (`genion`)

**입력 파일:** `ions.mdp`
```ini
integrator  = steep
nsteps      = 0
pbc         = xyz
```

```bash
$GMX grompp -f ions.mdp -c nanobody_solv.gro -p topol.top -o ions.tpr
$GMX genion -s ions.tpr -o nanobody_solv_ions.gro -p topol.top \
           -pname NA -nname CL -neutral
```

**출력:** `nanobody_solv_ions.gro`
- **CL⁻ 4개** 추가 → 전하 중성화 (0e)
- 총 잔기: 125 Protein + 7442 Water + 4 Ion

---

## Step 6: 에너지 최소화 (`mdrun`)

**입력 파일:** `em.mdp`
```ini
integrator  = steep       ; steepest descent
emtol       = 1000.0      ; 수렴 기준 (kJ/mol/nm)
emstep      = 0.01
nsteps      = 50000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
```

```bash
$GMX grompp -f em.mdp -c nanobody_solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -v -deffnm em -ntmpi 1 -ntomp 4
```

**출력:** `em.gro`, `em.edr`, `em.log`, `em.trr`
- 최종 Potential energy: **-352,202 kJ/mol** ✅

---

## Step 7: NVT Equilibration

**입력 파일:** `nvt.mdp`
```ini
define          = -DPOSRES        ; 단백질 위치 고정
integrator      = md
nsteps          = 50000           ; 100 ps
dt              = 0.002           ; 2 fs
constraints     = h-bonds
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
ref_t           = 300   300       ; 300 K
pcoupl          = no
gen_vel         = yes
gen_temp        = 300
```

```bash
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
$GMX mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu
```

**출력:** `nvt.gro`, `nvt.cpt`, `nvt.edr`, `nvt.log`
- **Performance: 73.942 ns/day** ✅ (완료)

---

## Step 8: NPT Equilibration

**입력 파일:** `npt.mdp`
```ini
define          = -DPOSRES
integrator      = md
nsteps          = 50000           ; 100 ps
continuation    = yes
tcoupl          = V-rescale
ref_t           = 300   300
pcoupl          = Parrinello-Rahman
tau_p           = 2.0
ref_p           = 1.0             ; 1 bar
compressibility = 4.5e-5
gen_vel         = no
```

```bash
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu
```

**출력:** `npt.gro`, `npt.cpt`, `npt.edr`, `npt.log`
- **Performance: 74.500 ns/day** ✅ (완료)

---

## Step 9: Production MD (100 ns)

**입력 파일:** `md.mdp`
```ini
integrator      = md
nsteps          = 50000000        ; 100 ns (2 fs × 50M steps)
dt              = 0.002
nstxout-compressed = 5000         ; 10 ps마다 좌표 저장
nstenergy       = 5000
continuation    = yes
constraints     = h-bonds
tcoupl          = V-rescale
ref_t           = 300   300
pcoupl          = Parrinello-Rahman
ref_p           = 1.0
gen_vel         = no
```

```bash
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2
$GMX mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu
```

**출력:** `md.xtc`, `md.edr`, `md.log`, `md.cpt`
- 저장 간격: 10 ps → 총 **10,000 frames**
- GPU PME 가속 사용 (`-pme gpu`)

**자동화 스크립트 (`run_md.sh`):**
```bash
bash run_md.sh > md_run.log 2>&1 &
echo $! > md.pid
```

---

## Step 10: WT gmx_MMPBSA (MD 완료 후)

**실행 환경:** `conda activate docking-md`

> **gmx_MMPBSA 제약사항:** `&alanine_scanning`은 ALA/GLY 치환만 지원하며, PRO 원잔기는 제외됨.
> → WT 절대 안정성 계산 후, 13개 mutant는 별도 10 ns MD를 실행하는 전략 채택

### 10-1. Trajectory 전처리
```bash
# PBC 보정 (시스템 GMX 사용 — md.tpr 생성 버전과 일치)
/usr/local/gromacs/bin/gmx trjconv \
    -s md.tpr -f md.xtc -o md_noPBC.xtc \
    -pbc mol -center <<< $'1\n0'
```

### 10-2. WT MMPBSA 입력 파일 (`mmpbsa_stability.in`)
```ini
Input file for gmx_MMPBSA - Nanobody Stability (WT)
&general
   sys_name = "Nanobody_AF3",
   startframe = 1,
   endframe = 2000,
   interval = 5,       ; 400 frames 사용
   temperature = 300,
   forcefields = "oldff/leaprc.ff99SB-ILDN",
/
&gb
   igb = 2,            ; AMBER GBSA GB 모델 2
   saltcon = 0.150,    ; 생리적 염 농도 (150 mM)
/
```

> **공정한 비교를 위한 설정:** 13개 mutant는 10 ns MD를 실행하므로,
> WT도 마지막 10 ns 구간만 사용: `startframe = 9001, endframe = 10000`

### 10-3. WT 실행
```bash
conda activate docking-md
cd /home/laugh/rosetta/nanobody_humanization/gromacs_mmpbsa
bash run_mmpbsa.sh
```

**출력:** `mmpbsa_WT/FINAL_RESULTS.dat` — WT 절대 안정성 ΔG (kcal/mol)

---

## Step 11: 13개 Mutant MD + MMPBSA

### 전략
- WT: 100 ns MD의 **마지막 10 ns** 사용
- Mutant 13개: 각각 **10 ns MD** 실행 후 MMPBSA
- ΔΔG = ΔG(mutant) - ΔG(WT) — 동일 trajectory 길이로 공정 비교

**예상 소요 시간:** ~3.5시간/개 × 13개 = **약 45시간** (순차 실행)

### 11-1. Mutant PDB 생성 (`make_mutants.py`)

**pdbfixer** 사용 — backbone 유지 + 새 sidechain 자동 재구성

```bash
conda run -n docking-md python3 make_mutants.py
```

**Mutation 목록 (pdbfixer 형식: `ORIGRES-RESNUM-NEWRES`):**

| 이름 | pdbfixer 문자열 | 순번 |
|------|----------------|------|
| Q1E  | `GLN-1-GLU`    | 1 |
| Q5V  | `GLN-5-VAL`    | 5 |
| S11L | `SER-11-LEU`   | 11 |
| A14P | `ALA-14-PRO`   | 14 |
| G35S | `GLY-35-SER`   | 35 |
| V40A | `VAL-40-ALA`   | 40 |
| A49S | `ALA-49-SER`   | 49 |
| I50V | `ILE-50-VAL`   | 50 |
| A74S | `ALA-74-SER`   | 74 |
| K86R | `LYS-86-ARG`   | 86 |
| P87A | `PRO-87-ALA`   | 87 |
| M92V | `MET-92-VAL`   | 92 |
| Q120L| `GLN-120-LEU`  | 120 |

**출력:** `mutants/<NAME>/nanobody_<NAME>.pdb`

### 11-2. 단일 Mutant 파이프라인 (`run_single_mutant.sh`)

각 mutant에 대해 전체 파이프라인 실행:

```
pdb2gmx → editconf → solvate → genion → EM → NVT(100ps) → NPT(100ps) → MD(10ns) → MMPBSA
```

```bash
bash run_single_mutant.sh Q1E   # 단독 실행
```

**10 ns MD 파라미터 (`md_10ns.mdp`):**
- `nsteps = 5,000,000` (2 fs × 5M = 10 ns)
- `nstxout-compressed = 5000` → 1,000 frames (10 ps 간격)

### 11-3. 전체 자동화 실행 (`run_all_mutants.sh`)

```bash
# WT MD 완료 + WT MMPBSA 완료 후 실행
bash run_all_mutants.sh > all_mutants.log 2>&1 &
```

실행 순서:
1. `make_mutants.py` — 13개 mutant PDB 생성
2. `run_single_mutant.sh` × 13 — 순차 실행
3. `summarize_all_mmpbsa.py` — ΔΔG 집계

### 11-4. 결과 요약 (`summarize_all_mmpbsa.py`)

```bash
python3 summarize_all_mmpbsa.py
```

**출력:** `mmpbsa_ddg_summary.csv`
```
ΔΔG = ΔG(mutant) - ΔG(WT)  [kcal/mol]
양수 → 불안정화, 음수 → 안정화
```

---

## 모니터링 명령어

```bash
# Production MD 현재 스텝 확인
grep -E "^\s+[0-9]+ " gromacs_mmpbsa/md.log | tail -3

# MD 속도 확인 (중반 이후)
grep "Performance" gromacs_mmpbsa/md.log | tail -3

# 파이프라인 전체 로그
tail -5 gromacs_mmpbsa/md_run.log

# 생성된 파일 목록 (최신순)
ls -lht gromacs_mmpbsa/*.gro gromacs_mmpbsa/*.cpt gromacs_mmpbsa/*.xtc 2>/dev/null

# 프로세스 확인
ps aux | grep "gmx mdrun" | grep -v grep
```

---

## 파일 구조

```
gromacs_mmpbsa/
├── nanobody_af3.pdb              ← 입력 구조 (CIF → PDB 변환)
├── nanobody_processed.gro        ← pdb2gmx 출력 (H 추가, SS결합)
├── nanobody_box.gro              ← 박스 생성
├── nanobody_solv.gro             ← 용매화
├── nanobody_solv_ions.gro        ← 이온 추가
├── topol.top                     ← 시스템 토폴로지
├── posre.itp                     ← Position restraints
├── ions.mdp                      ← genion 용 dummy mdp
├── em.mdp / em.gro               ← 에너지 최소화 ✅
├── nvt.mdp / nvt.gro             ← NVT equilibration ✅
├── npt.mdp / npt.gro             ← NPT equilibration ✅
├── md.mdp / md.xtc               ← Production MD (100 ns) 🔄
├── md_10ns.mdp                   ← Mutant용 10 ns MD 파라미터
├── run_md.sh                     ← WT MD 자동화 스크립트
├── mmpbsa_stability.in           ← WT MMPBSA 입력
├── run_mmpbsa.sh                 ← WT MMPBSA 실행 스크립트
├── make_mutants.py               ← pdbfixer로 13개 mutant PDB 생성
├── run_single_mutant.sh          ← 단일 mutant 전체 파이프라인
├── run_all_mutants.sh            ← 13개 mutant 일괄 실행
└── summarize_all_mmpbsa.py       ← WT + 13 mutant ΔΔG 비교 요약
```

**실행 후 생성되는 구조:**
```
gromacs_mmpbsa/
├── md_noPBC.xtc                  ← WT PBC 보정 trajectory
├── mmpbsa_WT/
│   ├── FINAL_RESULTS.dat         ← WT ΔG(stability)
│   └── FINAL_RESULTS.csv
├── mutants/
│   ├── Q1E/
│   │   ├── nanobody_Q1E.pdb      ← pdbfixer 생성 mutant PDB
│   │   ├── topol.top, posre.itp
│   │   ├── em.gro, nvt.gro, npt.gro
│   │   ├── md.xtc                ← 10 ns MD trajectory
│   │   ├── mmpbsa/
│   │   │   └── FINAL_RESULTS.dat ← mutant ΔG(stability)
│   │   └── pipeline.log
│   ├── Q5V/ ...
│   └── (13개 mutant 동일 구조)
├── mmpbsa_ddg_summary.csv        ← 최종 ΔΔG 비교표
└── all_mutants.log               ← 전체 파이프라인 로그
```

---

## 진행 현황

| 단계 | 상태 | 비고 |
|------|------|------|
| Step 1-5 (준비) | ✅ 완료 | |
| Step 6 (EM) | ✅ 완료 | -352,202 kJ/mol |
| Step 7 (NVT) | ✅ 완료 | 73.942 ns/day |
| Step 8 (NPT) | ✅ 완료 | 74.500 ns/day |
| Step 9 (MD 100 ns) | 🔄 실행 중 | PID: 3024538, ~75 ns/day, 완료까지 ~31시간 |
| Step 10 (WT MMPBSA) | ⏳ 대기 | `bash run_mmpbsa.sh` |
| Step 11 (Mutant 13개 MD+MMPBSA) | ⏳ 대기 | `bash run_all_mutants.sh`, ~45시간 |

---

**작성일:** 2026-02-10
**Force field:** AMBER99SB-ILDN
**Water model:** TIP3P
**MD 길이:** 100 ns
**분석 도구:** gmx_MMPBSA 1.6.4 (docking-md conda env)
**GROMACS:** 2024.1 (`/usr/local/gromacs/bin/gmx`, CUDA GPU)
