import pandas as pd
import os
import argparse

def read_data(file_path, sheet_name, save_path):
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    columns_to_drop = ['Mouse_ID', 'Tumor_Side', 'Age']
    df = df.drop(columns=columns_to_drop)
    treatment_type = [0, 1, 11, 10]
    sham_df = df[df['Treatment'] == treatment_type[0]].drop(columns='Treatment')
    ir_df = df[df['Treatment'] == treatment_type[1]].drop(columns='Treatment')
    aspirin_df = df[df['Treatment'] == treatment_type[2]].drop(columns='Treatment')
    ir_aspirin_df = df[df['Treatment'] == treatment_type[3]].drop(columns='Treatment')
    columns_path = os.path.join(save_path, "columns.txt")
    with open(columns_path, 'w') as f:
        f.write(','.join(df.columns))
    return sham_df, ir_df, aspirin_df, ir_aspirin_df

def combine_dfs(sham_df, ir_df, aspirin_df, ir_aspirin_df, save_path):
    dfs = {'sham': sham_df, 'ir': ir_df, 'aspirin': aspirin_df, 'ir+aspirin': ir_aspirin_df}
    merged_dataframes = []
    for name1, df1 in dfs.items():
        for name2, df2 in dfs.items():
            if name1 < name2:
                df1_modified = df1.copy()
                df1_modified['Origin'] = 0
                df1_modified = df1_modified[['Origin'] + [col for col in df1_modified.columns if col != 'Origin']]
                df2_modified = df2.copy()
                df2_modified['Origin'] = 1
                df2_modified = df2_modified[['Origin'] + [col for col in df2_modified.columns if col != 'Origin']]
                merged_df = pd.concat([df1_modified, df2_modified])
                merged_dataframes.append((f"{name1}_{name2}", merged_df))
    for name, df in merged_dataframes:
        full_path = os.path.join(save_path, f"{name}.csv")
        df.to_csv(full_path, index=False, header=False)
        print(f"Saved DataFrame '{name}' to '{full_path}'")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument('sheet_name', type=str, help='Excel sheet name to be processed')
    args = parser.parse_args()
    
    file_path = '.\\raw_data\\Final_Spreadsheet_separated_by_age.xlsx'
    save_path = '.\\data\\output_preprocess'
    sham_df, ir_df, aspirin_df, ir_aspirin_df = read_data(file_path, args.sheet_name, save_path)
    combine_dfs(sham_df, ir_df, aspirin_df, ir_aspirin_df, save_path)


## how to run the code:  python combine_dfs.py Age_W5

