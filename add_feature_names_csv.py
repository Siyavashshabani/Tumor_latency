import pandas as pd
import os

def add_feature_names_csv(names_file_path, csv_directory_path, save_directory_path):
    # Read names from the file into a list
    with open(names_file_path, 'r') as file:
        names = [line.strip() for line in file]

    # Read CSV files from directory
    csv_files = [f for f in os.listdir(csv_directory_path) if f.endswith('.csv')]

    # Function to map indices to names
    def index_to_name(index):
        # Check if index is an integer and within the valid range
        if isinstance(index, int) and 1 <= index <= len(names):
            return names[index - 1]  # Convert 1-based index to 0-based
        else:
            return 'Index out of range'

    # Process each CSV file
    for csv_file in csv_files:
        # Read the current CSV file
        current_table = pd.read_csv(os.path.join(csv_directory_path, csv_file))
        # Convert 'V1' to integers, handling non-convertible values by coercing to NaN
        current_table['V1'] = pd.to_numeric(current_table['V1'], errors='coerce').fillna(0).astype(int)
        # Map the 'V1' column to names using the provided function
        current_table['Names'] = current_table['V1'].apply(index_to_name)

        # Save the modified table
        modified_file_path = os.path.join(save_directory_path, 'modified_' + csv_file)
        current_table.to_csv(modified_file_path, index=False)
        print(f"Processed and saved modified data to {modified_file_path}")

if __name__ == "__main__":
    # Define paths
    names_file_path = 'C:\\Shabani\\Projects\\tumor_latency\\data\\column_names.csv'
    csv_directory_path = 'C:\\Shabani\\Projects\\tumor_latency\\data\\output_r\\'
    save_directory_path = 'C:\\Shabani\\Projects\\tumor_latency\\data\\output_r_joint_features\\'

    # Call the function
    add_feature_names_csv(names_file_path, csv_directory_path, save_directory_path)



## How to run the code: python add_feature_names_csv.py