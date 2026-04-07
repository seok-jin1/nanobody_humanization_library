#!/usr/bin/env python3
"""
숫자 IMGT 포지션 기반으로 나노바디 잔기 다양성(엔트로피) 및 아미노산 물성 비율을 시각화한다.

Visualize nanobody residue variability by numeric IMGT position (1-128):
upper panel shows Shannon entropy per position with CDR highlights,
lower panel shows stacked amino-acid property ratios.

Input:
    - 01_imgt_residue_matrix.csv : residue frequency ratio matrix
Output:
    - figures/04_imgt_numeric_only_analysis.png : entropy + property ratio plot

Usage:
    python 04_visualize_numeric_only.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- 설정 ---
INPUT_CSV = '01_imgt_residue_matrix.csv'
OUTPUT_PLOT = 'figures/04_imgt_numeric_only_analysis.png'

PROPERTIES = {
    'Hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'],
    'Polar': ['S', 'T', 'N', 'Q', 'Y'],
    'Positive': ['K', 'R', 'H'],
    'Negative': ['D', 'E'],
    'Special': ['C', 'G'],
    'Others': ['X', '-']
}

COLOR_MAP = {
    'Hydrophobic': '#e41a1c', 
    'Polar': '#377eb8',       
    'Positive': '#4daf4a',    
    'Negative': '#984ea3',    
    'Special': '#ff7f00',     
    'Others': '#999999'       
}

def calculate_entropy(col_series):
    """아미노산 분포에 대한 Shannon Entropy 계산"""
    probs = col_series[col_series > 0]
    if len(probs) <= 1: return 0
    # 합이 1이 되도록 정규화 (집계 과정에서 합이 1이 아닐 수 있음)
    probs = probs / probs.sum()
    return -np.sum(probs * np.log2(probs))

def get_num(pos):
    """문자열 포지션에서 숫자만 추출"""
    num_str = "".join([c for c in pos if c.isdigit()])
    return int(num_str) if num_str else 0

def main():
    # 1. 데이터 로드
    df = pd.read_csv(INPUT_CSV, index_col=0)
    df = df.groupby(df.index).sum() # 인덱스 중복 방지
    
    # 2. 숫자 부분으로 그룹화하여 집계 (Mean 사용)
    # 컬럼명을 숫자로 변환하여 그룹화
    df.columns = [get_num(c) for c in df.columns]
    
    # 같은 숫자를 가진 컬럼들을 평균(Mean) 내어 하나의 대표 프로파일로 만듦
    # (주의: Sum을 하면 비율이 1을 넘을 수 있으므로 Mean이 적절함)
    df_numeric = df.groupby(level=0, axis=1).mean()
    
    # 1-128 범위로 필터링 및 정렬
    numeric_positions = sorted([c for c in df_numeric.columns if 1 <= c <= 128])
    df_numeric = df_numeric[numeric_positions]
    
    # 3. 데이터 계산
    entropy_list = []
    prop_data = {k: [] for k in PROPERTIES.keys()}
    
    for pos in numeric_positions:
        col = df_numeric[pos]
        entropy_list.append(calculate_entropy(col))
        for prop, aas in PROPERTIES.items():
            available_aas = [aa for aa in aas if aa in col.index]
            prop_sum = col.loc[available_aas].sum()
            prop_data[prop].append(prop_sum)

    # 4. 시각화
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 11), sharex=True)

    # 상단: Entropy Plot
    ax1.plot(numeric_positions, entropy_list, color='black', linewidth=2, label='Shannon Entropy')
    ax1.fill_between(numeric_positions, entropy_list, color='black', alpha=0.1)
    ax1.set_ylabel('Entropy (bits)', fontsize=12)
    ax1.set_title('Nanobody Variability by Numeric IMGT Position (1-128)', fontsize=18)

    # CDR 하이라이트 (숫자 기준이므로 더 정확함)
    def highlight_cdr(ax):
        cdrs = [(27, 38, 'CDR1'), (56, 65, 'CDR2'), (105, 117, 'CDR3')]
        colors = ['#ffcccc', '#ccffcc', '#ccccff']
        for (start, end, label), color in zip(cdrs, colors):
            ax.axvspan(start, end, color=color, alpha=0.4, label=label if ax == ax1 else "")
            if ax == ax1:
                ax.text((start+end)/2, max(entropy_list)*1.05, label, ha='center', fontweight='bold')

    highlight_cdr(ax1)
    highlight_cdr(ax2)

    # 하단: Stacked Bar Chart
    bottom = np.zeros(len(numeric_positions))
    for prop, values in prop_data.items():
        # 비율 정규화 (바의 높이를 1로 맞춤)
        vals = np.array(values)
        ax2.bar(numeric_positions, vals, bottom=bottom, label=prop, color=COLOR_MAP[prop], width=0.8)
        bottom += vals

    ax2.set_ylabel('Property Ratio', fontsize=12)
    ax2.set_xlabel('IMGT Numeric Position', fontsize=12)
    
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right', bbox_to_anchor=(1.08, 1.0))

    # X축 눈금을 5단위로 표시하여 깔끔하게 만듦
    plt.xticks(np.arange(0, 131, 5))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Numeric-only visualization saved to {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
