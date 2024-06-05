import pandas as pd
import os
import argparse

def main(age):
    # Set the directory path for the CSV files
    csv_directory_path = 'C:\\Shabani\\Projects\\tumor_latency\\data\\output_r_joint_features\\'

    # Set the path for the names
    names_path = "C:\\Shabani\\Projects\\tumor_latency\\data\\names.csv"
    
    # Set the directory path where the results should be saved
    saved_path = 'C:\\Shabani\\Projects\\tumor_latency\\data\\analysis\\'

    # Read the names from the column table
    names = pd.read_csv(names_path)

    # List all CSV files in the directory
    csv_files = [f for f in os.listdir(csv_directory_path) if f.endswith('.csv')]

    # Initialize the result DataFrame with names
    result_df = names.copy()

    # Loop through each file in the csv_files list
    for csv_file in csv_files:
        # Read the current p-value table
        current_table = pd.read_csv(os.path.join(csv_directory_path, csv_file))
        
        # Extract the base name without the '.csv' extension for use in naming the column
        table_name = csv_file.replace('.csv', '')
        
        # Merge the names DataFrame with the current table on the 'Names' column
        temp_df = pd.merge(names[['Names']], current_table[['Names', 'V7']], on='Names', how='left')
        
        # Rename the 'V7' column to include the name of the current CSV file, slicing from index 16 if needed
        temp_df = temp_df.rename(columns={'V7': f'{table_name[16:]}'})

        # Add the new column to the result DataFrame
        result_df[f'{table_name[16:]}'] = temp_df[f'{table_name[16:]}']

    # Specify the file path where the result should be saved
    save_result_path = os.path.join(saved_path, f'analysis_{age}.csv')

    # Save the result DataFrame to a CSV file
    result_df.to_csv(save_result_path, index=False)

    # Print a confirmation message
    print(f"Results have been successfully saved to {save_result_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process CSV files and merge data based on column names.")
    parser.add_argument('age', type=str, help="Age for the output file name to identify the analysis result.")
    args = parser.parse_args()
    
    main(args.age)
