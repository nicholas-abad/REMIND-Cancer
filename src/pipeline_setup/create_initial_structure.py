import os
import shutil
import glob
import json
import pandas as pd
import argparse
from tqdm import tqdm

def create_folder_structure(metadata_path, output_folder):
    """
    Creates a structured folder system for the pipeline.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata CSV file.
    output_folder : str
        Output path where patient subfolders will be created.
    """
    # Load metadata
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    metadata = pd.read_csv(metadata_path)

    # Create main output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    num_vcf_not_copied = 0

    for _, row in tqdm(metadata.iterrows(), total=len(metadata), desc="Setting up folder structure"):
        pid, tumor_origin, vcf_path = row["pid"], row["tumor_origin"], row["path_to_wgs_file"]

        if not os.path.exists(vcf_path):
            num_vcf_not_copied += 1
            print(f"Missing VCF: {vcf_path}")
            continue

        patient_folder = os.path.join(output_folder, f"{pid}_{tumor_origin}")
        os.makedirs(patient_folder, exist_ok=True)

        new_vcf_name = os.path.join(patient_folder, os.path.basename(vcf_path).replace(".vcf", "_original.vcf"))
        if not os.path.exists(new_vcf_name):
            shutil.copy(vcf_path, new_vcf_name)

    print(f"\n✅ Folder structure created at: {output_folder}")
    print(f"❗ Number of VCF files not copied: {num_vcf_not_copied}")


def create_results_json(metadata_path, patient_folders_path, results_json_path):
    """
    Generates a results tracking JSON file.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata CSV file.
    patient_folders_path : str
        Path to the folder containing all patient subfolders.
    results_json_path : str
        Output path for the results JSON file.
    """
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    metadata = pd.read_csv(metadata_path)

    results_dict = {
        "results": {
            "original": {
                "primary_tumor_wgs": [],
                "primary_tumor_wgs_and_rnaseq": [],
                "metastasic_tumor_wgs": [],
                "metastasic_tumor_wgs_and_rnaseq": []
            }
        }
    }

    for folder in tqdm(glob.glob(os.path.join(patient_folders_path, "*")), desc="Generating JSON file"):
        vcf_files = glob.glob(os.path.join(folder, "*_original.vcf"))
        if not vcf_files:
            continue

        vcf_file = vcf_files[0]
        pid = os.path.basename(folder).split("_")[0]

        if pid not in metadata["pid"].values:
            raise ValueError(f"{pid} does not exist in the metadata dataframe.")

        tumor_origin = metadata.loc[metadata["pid"] == pid, "tumor_origin"].values[0]

        ge_data_available = False
        if "path_to_rnaseq_file" in metadata.columns:
            ge_data_available = os.path.exists(metadata.loc[metadata["pid"] == pid, "path_to_rnaseq_file"].values[0])
        elif "ge_data_available" in metadata.columns:
            ge_data_available = metadata.loc[metadata["pid"] == pid, "ge_data_available"].values[0]

        if tumor_origin == "primary_tumor":
            category = "primary_tumor_wgs_and_rnaseq" if ge_data_available else "primary_tumor_wgs"
        else:
            category = "metastasic_tumor_wgs_and_rnaseq" if ge_data_available else "metastasic_tumor_wgs"

        results_dict["results"]["original"][category].append(vcf_file)

    with open(results_json_path, "w") as f:
        json.dump(results_dict, f, indent=4)

    print(f"\n✅ Results JSON saved at: {results_json_path}")


def main(config_path):
    """
    Main function to set up folder structure and create the results tracking JSON.

    Parameters
    ----------
    config_path : str
        Path to the configuration JSON file.
    """
    # Load configuration file
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, "r") as f:
        config = json.load(f)

    dataset = config["pipeline"]["dataset"]
    metadata_path = config["pipeline"]["path_to_metadata"]
    patient_folders_path = config["pipeline_folder_setup"][dataset]["output_path_to_patient_folders"]
    results_json_path = config["pipeline"]["path_to_results"]

    print("\n#### (1/2) Creating Folder Structure ####")
    create_folder_structure(metadata_path, patient_folders_path)

    print("\n#### (2/2) Creating Results JSON File ####")
    create_results_json(metadata_path, patient_folders_path, results_json_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setup REMIND-Cancer folder structure and create results JSON file.")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the configuration JSON file.")

    args = parser.parse_args()
    main(args.config)