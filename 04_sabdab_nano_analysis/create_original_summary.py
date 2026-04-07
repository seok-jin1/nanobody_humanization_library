import pandas as pd

# --- 설정 ---
INPUT_CSV = 'imgt_residue_matrix.csv'
OUTPUT_SUMMARY_CSV = 'imgt_original_summary_top20.csv'

def main():
    # 1. 데이터 로드
    # index_col=0 은 아미노산 잔기들 (A, R, N, D, ..., -, X)
    df = pd.read_csv(INPUT_CSV, index_col=0)
    
    # 인덱스 중복 방지 (만약 생성 과정에서 발생했을 경우 대비)
    df = df.groupby(df.index).sum()
    
    # 2. 모든 포지션(원본 컬럼명 그대로) 순회
    summary_data = []
    
    for pos in df.columns:
        col = df[pos]
        # 20% 이상인 것 필터링
        top_residues = col[col >= 0.2].sort_values(ascending=False)
        
        if top_residues.empty:
            continue
            
        # 문자열 생성: "AA(Ratio), AA(Ratio)"
        res_str = ", ".join([f"{aa}({ratio:.2f})" for aa, ratio in top_residues.items()])
        
        summary_data.append({
            'IMGT_Position': pos,
            'Top_Residues_20pct': res_str
        })
    
    # 3. CSV 저장
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(OUTPUT_SUMMARY_CSV, index=False)
    print(f"Original position summary CSV (top 20%) saved to {OUTPUT_SUMMARY_CSV}")

if __name__ == "__main__":
    main()
