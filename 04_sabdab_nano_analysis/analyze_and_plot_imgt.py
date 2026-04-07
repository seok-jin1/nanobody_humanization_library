import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- 설정 ---
SUMMARY_FILE = 'sabdab_nano_summary_all.tsv'
IMGT_DIR = 'imgt'
OUTPUT_CSV = 'imgt_residue_matrix.csv'
OUTPUT_PLOT = 'imgt_residue_heatmap.png'

# 3-letter to 1-letter mapping
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X' 
}

def load_chain_mapping(summary_file):
    mapping = {}
    try:
        df = pd.read_csv(summary_file, sep='\t')
        for _, row in df.iterrows():
            pdb = str(row['pdb']).lower()
            hchain = str(row['Hchain'])
            if hchain and hchain != 'nan':
                mapping[pdb] = hchain
    except Exception as e:
        print(f"Error reading summary file: {e}")
    return mapping

def parse_pdb_residues(pdb_path, target_chain):
    residues = {} 
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    if chain_id != target_chain:
                        continue
                    atom_name = line[12:16].strip()
                    if atom_name != 'CA':
                        continue
                    res_name = line[17:20].strip()
                    res_seq = line[22:26].strip()
                    i_code = line[26].strip()
                    # IMGT format: position + optional insertion code
                    imgt_pos = res_seq + i_code
                    aa = AA_MAP.get(res_name, 'X')
                    residues[imgt_pos] = aa
    except Exception as e:
        print(f"Error parsing {pdb_path}: {e}")
        return None
    return residues

def sort_imgt_positions(positions):
    unique_positions = sorted(list(set(positions)))
    def key_func(pos):
        num_part = "".join([c for c in pos if c.isdigit()])
        char_part = "".join([c for c in pos if not c.isdigit()])
        return (int(num_part) if num_part else 0, char_part)
    return sorted(unique_positions, key=key_func)

def main():
    print("1. Loading metadata...")
    chain_map = load_chain_mapping(SUMMARY_FILE)

    pdb_files = glob.glob(os.path.join(IMGT_DIR, '*.pdb'))
    print(f"2. Processing {len(pdb_files)} PDB files...")

    data = []
    valid_pdbs = 0
    for pdb_path in pdb_files:
        pdb_id = os.path.basename(pdb_path)[:4].lower()
        if pdb_id not in chain_map:
            continue
        target_chain = chain_map[pdb_id]
        residues = parse_pdb_residues(pdb_path, target_chain)
        if residues:
            data.append(residues)
            valid_pdbs += 1

    if valid_pdbs == 0:
        print("No valid data found.")
        return

    # DataFrame 생성 (각 행은 하나의 PDB)
    df = pd.DataFrame(data)
    df.fillna('-', inplace=True)
    
    # 중복 컬럼 제거 (혹시 있다면)
    df = df.loc[:, ~df.columns.duplicated()]
    
    # IMGT Position 정렬 및 필터링
    all_cols = df.columns.tolist()
    sorted_cols = sort_imgt_positions(all_cols)
    
    print("3. Calculating frequency matrix...")
    aa_list = sorted(list(set(AA_MAP.values()))) + ['X', '-']
    # Matrix를 먼저 빈 값으로 만들고 채움
    matrix = pd.DataFrame(0, index=aa_list, columns=sorted_cols)
    
    for col in sorted_cols:
        if col in df.columns:
            counts = df[col].value_counts()
            for aa, count in counts.items():
                if aa in matrix.index:
                    matrix.at[aa, col] = count

    ratio_matrix = matrix / valid_pdbs
    ratio_matrix.to_csv(OUTPUT_CSV)

    print("4. Generating visualization...")
    def get_num(pos):
        num_str = "".join([c for c in pos if c.isdigit()])
        return int(num_str) if num_str else 0

    # Nanobody V-domain (1-128)
    v_domain_cols = [col for col in sorted_cols if 1 <= get_num(col) <= 128]
    
    # Y축 순서 (Gap을 맨 아래로)
    y_order = sorted([aa for aa in aa_list if aa not in ['X', '-']]) + ['X', '-']
    
    # Reindex 전에 matrix 중복 여부 확인
    plot_data = ratio_matrix.loc[y_order, v_domain_cols]

    plt.figure(figsize=(22, 8))
    sns.heatmap(plot_data, cmap='YlGnBu', cbar_kws={'label': 'Frequency Ratio'})
    
    plt.title(f'Nanobody Residue Distribution (IMGT 1-128, n={valid_pdbs})')
    plt.xlabel('IMGT Position')
    plt.ylabel('Amino Acid')
    
    # X축 눈금이 너무 많으면 읽기 힘드므로 적절히 조절
    if len(v_domain_cols) > 100:
        plt.xticks(rotation=90, fontsize=6)
    else:
        plt.xticks(rotation=90, fontsize=8)
        
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"   Saved heatmap to {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
