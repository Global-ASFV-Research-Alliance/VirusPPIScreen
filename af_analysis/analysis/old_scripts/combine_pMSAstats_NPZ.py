import pandas as pd
import sys
import os

def combine_dataframes_columns(file_path1, file_path2):
    # Load the dataframes from the CSV files
    df1 = pd.read_csv(file_path1, index_col=0)
    df2 = pd.read_csv(file_path2, index_col=0)

    # Perform an outer join to combine the columns from both dataframes
    combined_df = pd.merge(df1, df2, left_index=True, right_index=True, how='outer')

    return combined_df

if __name__ == "__main__":
    # Take file paths from command-line arguments
    if len(sys.argv) != 4:
        raise ValueError("Please provide two file paths for the dataframes.")

    input_file_path1 = sys.argv[1]
    input_file_path2 = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Combine the columns from the two dataframes
    combined_df = combine_dataframes_columns(input_file_path1, input_file_path2)

    # Create the output file path with a '_combined' tag
    basename = os.path.basename(input_file_path1)
    output_file_path = f"{output_dir}/{os.path.splitext(basename)[0]}_combined.csv"

    # Save the combined dataframe to a CSV file
    combined_df.to_csv(output_file_path)
    print(f"Combined dataframe saved to {output_file_path}")
