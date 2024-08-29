import os
import shutil
import glob
import pandas as pd
from tqdm import tqdm as tqdm
import json

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import _get_pid_from_structured_vcf_path


def create_new_structure(
    path_to_metadata_dataframe: str,
    path_to_patient_folders: str,
    ending_of_original_files: str,    
):
    """Create a folder structure that will fit with the pipeline.

    Parameters
    ----------
    path_to_metadata_dataframe : str
        Path to the metadata dataframe that was previously created in ./PCAWG/*.
    path_to_patient_folders : str
        Output path to the folder that will contain subfolders for each patient.
        Example of Final Structure (MASTER):
            path_to_patient_folders
                pcawg_pid_1
                    snvs_pcawg_pid1_somatic_snvs_conf_8_to_10.vcf
                pcawg_pid_2
                    snvs_pcawg_pid2_somatic_snvs_conf_8_to_10.vcf
    """
    print("Creating new structure...")
    # Read in the previously-created metadata dataframe.
    data = pd.read_csv(path_to_metadata_dataframe)
 
    # Create the initial folder (path_to_patient_folders) if not yet created. 
    if not os.path.exists(path_to_patient_folders):
        os.mkdir(path_to_patient_folders)

    # For each patient in the metadata dataframe, create a subfolder and copy the original SNV .vcf into the patient's subfolder.
    for _, row in tqdm(data.iterrows(), total=data.shape[0], desc="Copying SNV Files"):
        pid = row["pid"]
        cancer_type = row["cancer_type"]
        path_to_vcf = row["path_to_vcf"]
        
        if not os.path.exists(path_to_vcf):
            print(f"{path_to_vcf} does not exist. Skipping...")
            continue
        
        if (path_to_vcf != "not_available") and (os.path.exists(path_to_vcf)):
            # Create folder for single patient if this does not exist.
            path_to_single_patient_folder = os.path.join(path_to_patient_folders, f"{pid}_{cancer_type}")
            if not os.path.exists(path_to_single_patient_folder):
                os.mkdir(path_to_single_patient_folder)
                
            # Save .vcf file to the folder for the single patient if this does not exist.
            if not os.path.exists(os.path.join(path_to_single_patient_folder, path_to_vcf.split("/")[-1])):
                name_of_original_vcf = path_to_vcf.split("/")[-1]
                name_of_new_vcf = os.path.join(path_to_single_patient_folder, name_of_original_vcf.replace(".vcf", ending_of_original_files))
                if not os.path.exists(name_of_new_vcf):
                    shutil.copy(path_to_vcf, name_of_new_vcf)

    print(f".... Folder structure completed: {path_to_patient_folders} ")

def create_json_file(
        path_to_metadata_dataframe: str,
        path_to_patient_folders: str,
        path_to_results_json_file: str,
        ending_of_original_files: str,    
):
    results_dict = {
        "results": {
            "original": {
                "tumor_snv": [],
                "tumor_snv_and_ge": [],
                "metastasis_snv": [],
                "metastasis_snv_and_ge": []
            }
        }
    }
    metadata = pd.read_csv(path_to_metadata_dataframe)
    for folder in tqdm(glob.glob(os.path.join(path_to_patient_folders, "*"))):
        # Get the original .vcf filename.
        original_file = glob.glob(os.path.join(folder, f"*{ending_of_original_files}"))  # Returns list.
        original_file = original_file[0]
        
        # Get the PID.
        pid = _get_pid_from_structured_vcf_path(original_file, only_pid=True)
        assert pid in list(metadata["pid"]), f"{pid} does not exist in metadata dataframe."
        
        # Get the cancer_type (tumor or metastasis) and ge_data_available (True or False).
        cancer_type = metadata[metadata["pid"] == pid].iloc[0]["cancer_type"]

        if "path_to_tsv" in list(metadata.columns):
            ge_data_available = True if metadata[metadata["pid"] == pid].iloc[0]["path_to_tsv"] != "not_available" else False
        else:
            ge_data_available = metadata[metadata["pid"] == pid].iloc[0]["ge_data_available"]
        
        if (cancer_type == "tumor"):
            if (ge_data_available):
                results_dict["results"]["original"]["tumor_snv_and_ge"].append(original_file)
            else:
                results_dict["results"]["original"]["tumor_snv"].append(original_file)
        else:
            if (ge_data_available):
                results_dict["results"]["original"]["metastasis_snv_and_ge"].append(original_file)
            else:
                results_dict["results"]["original"]["metastasis_snv"].append(original_file)
                
    with open(path_to_results_json_file, "w") as f:
        json.dump(results_dict, f)
        
        