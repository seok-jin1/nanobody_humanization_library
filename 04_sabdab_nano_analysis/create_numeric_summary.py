import pandas as pd
import numpy as np

# --- 설정 ---
INPUT_CSV = 'imgt_residue_matrix.csv'
OUTPUT_SUMMARY_CSV = 'imgt_numeric_summary_top20.csv'

def get_num(pos):
    """문자열 포지션에서 숫자만 추출"""
    num_str = "".join([c for c in pos if c.isdigit()])
    return int(num_str) if num_str else 0

def main():
    # 1. 데이터 로드
    df = pd.read_csv(INPUT_CSV, index_col=0)
    
    # 2. 숫자 부분으로 그룹화하여 평균 계산
    df.columns = [get_num(c) for c in df.columns]
    df_numeric = df.groupby(level=0, axis=1).mean()
    
    # 1-128 범위 제한
    numeric_positions = sorted([c for c in df_numeric.columns if 1 <= c <= 128])
    df_numeric = df_numeric[numeric_positions]
    
    # 3. 20% 이상인 잔기 추출 및 포맷팅
    summary_data = []
    
    for pos in numeric_positions:
        col = df_numeric[pos]
        # 20% 이상인 것 필터링 (Gap '-' 포함 여부는 데이터 성격에 따라 결정되나 여기선 포함)
        top_residues = col[col >= 0.2].sort_values(ascending=False)
        
        # 문자열 생성: "AA(Ratio), AA(Ratio)"
        res_str = ", ".join([f"{aa}({ratio:.2f})" for aa, ratio in top_residues.items()])
        
        summary_data.append({
            'IMGT_Position': pos,
            'Top_Residues_20pct': res_str
        })
    
    # 4. CSV 저장
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(OUTPUT_SUMMARY_CSV, index=False)
    print(f"Summary CSV (top 20%) saved to {OUTPUT_SUMMARY_CSV}")

if __name__ == "__main__":
    main()
