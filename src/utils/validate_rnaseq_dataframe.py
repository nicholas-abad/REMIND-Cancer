import os
import pandas as pd
import argparse

def validate_rnaseq_dataframe(rnaseq_path):
    """
    Validates an RNA-seq dataframe by checking the presence of 'pid' column, ensuring the delimiter is a comma,
    and counting the number of PIDs and genes.

    Parameters
    ----------
    rnaseq_path : str
        Path to the RNA-seq dataframe CSV file.
    """
    # Ensure the RNA-seq file exists
    if not os.path.exists(rnaseq_path):
        print(f"❌ RNA-seq file not found: {rnaseq_path}")
        return

    # Load RNA-seq dataframe
    print("\n🔍 (1/2) Attempting to load in RNA-Seq dataframe...")
    try:
        df = pd.read_csv(rnaseq_path, sep=",")
        print(f"✅ RNA-Seq dataframe can be loaded in properly.")
    except Exception as e:
        print(f"❌ Error reading RNA-seq file: {e}")
        return

    # Check if 'pid' column exists
    print("\n🔍 (2/2) Checking if 'pid' is a column in the dataframe...")
    if "pid" not in df.columns:
        print("❌ RNA-seq file is missing the required 'pid' column.")
        return

    # Count PIDs and genes
    num_pids = len(df)
    num_genes = len(df.columns) - 1  # Exclude 'pid' column

    print(f"✅ RNA-seq file is valid: {num_pids} PIDs and {num_genes} genes detected.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate RNA-seq dataframe format.")
    parser.add_argument("-r", "--rnaseq", required=True, type=str, help="Path to the RNA-seq dataframe CSV file.")
    
    args = parser.parse_args()
    
    validate_rnaseq_dataframe(args.rnaseq)