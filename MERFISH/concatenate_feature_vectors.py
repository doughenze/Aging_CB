import pandas as pd
import glob
import argparse

def concatenate_files(output_dir, final_output):
    csv_files = glob.glob(f"{output_dir}/feature_vectors_*.csv")
    
    dfs = [pd.read_csv(f) for f in csv_files]
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(final_output, index=False)
    print(f"Concatenated file saved to {final_output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenate feature vector CSV files.')
    parser.add_argument('output_dir', type=str, help='Directory containing the output CSV files')
    parser.add_argument('final_output', type=str, help='Path to the final concatenated CSV file')
    args = parser.parse_args()
    
    concatenate_files(args.output_dir, args.final_output)
