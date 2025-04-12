import os
import pandas as pd
import argparse

def validate_vcf(vcf_path):
    """
    Validates a VCF file by checking its required columns.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.
    """
    required_columns = {"#CHROM", "POS", "REF", "ALT", "SEQUENCE_CONTEXT", "GENE"}
    
    # Ensure the VCF file exists
    if not os.path.exists(vcf_path):
        return False, f"❌ VCF file not found: {vcf_path}"
    
    # Try reading the VCF file
    try:
        df = pd.read_csv(vcf_path, sep="\t")
    except Exception as e:
        return False, f"❌ Error reading VCF file {vcf_path}: {e}"
    
    # Check if all required columns exist
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        return False, f"🚨 Missing columns in {vcf_path}: {', '.join(missing_columns)}"
    
    return True, None

def validate_metadata(metadata_path):
    """
    Validates the metadata CSV file by checking its required columns and verifying file existence.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata CSV file.
    """
    required_columns = {"pid", "tumor_origin", "path_to_wgs_file", "cohort", "ge_data_available"}

    # Validate metadata file.
    print("\n🔍 (1/2) Validating metadata file...")
    
    # Ensure the metadata file exists
    if not os.path.exists(metadata_path):
        print(f"❌ Metadata file not found: {metadata_path}")
        return
    else:
        print(f"✅ Metadata file exists.")

    # Load metadata CSV
    try:
        metadata_df = pd.read_csv(metadata_path)
        print(f"✅ Metadata can be loaded in properly.")
    except Exception as e:
        print(f"❌ Error reading metadata file: {e}")
        return

    # Check if all required columns exist
    missing_columns = required_columns - set(metadata_df.columns)
    if missing_columns:
        print(f"❌ Metadata file is missing required columns: {', '.join(missing_columns)}")
        return
    else:
        print(f"✅ Metadata contains all necessary columns.")

    # Check if WGS files exist and validate associated VCF files
    total_wgs_files = len(metadata_df)
    existing_wgs_files = metadata_df[metadata_df['path_to_wgs_file'].apply(lambda x: os.path.exists(x))]
    missing_wgs_files = metadata_df[~metadata_df['path_to_wgs_file'].apply(lambda x: os.path.exists(x))]
    
    num_existing = len(existing_wgs_files)
    num_missing = len(missing_wgs_files)

    print(f"✅ {num_existing} of {total_wgs_files} WGS files exist.")
    if num_missing > 0:
        print(f"🚨 {num_missing} WGS files are missing. They are:")
        for path in missing_wgs_files['path_to_wgs_file']:
            print(f"   - {path}")
    
    # Validate associated VCF files
    print("\n🔍 (2/2) Validating associated VCF files...")
    total_vcf_files = len(existing_wgs_files)
    valid_vcf_count = 0
    invalid_vcf_messages = []
    
    for wgs_file in existing_wgs_files['path_to_wgs_file']:
        vcf_path = wgs_file.replace(".wgs", ".vcf")  # Assuming VCF file has the same name but with .vcf extension
        is_valid, error_message = validate_vcf(vcf_path)
        if is_valid:
            valid_vcf_count += 1
        else:
            invalid_vcf_messages.append(error_message)
    
    print(f"✅ {valid_vcf_count} of {total_vcf_files} VCF files are valid.")
    if invalid_vcf_messages:
        print("\n🚨 Invalid VCF files detected:")
        for msg in invalid_vcf_messages:
            print(msg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate metadata CSV file format and WGS file existence.")
    parser.add_argument("-m", "--metadata", required=True, type=str, help="Path to the metadata CSV file.")
    
    args = parser.parse_args()
    
    validate_metadata(args.metadata)