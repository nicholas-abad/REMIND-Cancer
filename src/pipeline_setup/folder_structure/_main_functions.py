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
    output_path_to_patient_folders: str,
):
    """Create a folder structure that will fit with the pipeline.

    Parameters
    ----------
    path_to_metadata_dataframe : str
        Path to the previously-created metadata dataframe.
    output_path_to_patient_folders : str
        Output path to the folder that will contain subfolders for each patient.
        Example of Final Structure:
            output_path_to_patient_folders
                pid_1
                    wgs_snv_file_pid1_original.vcf
                pid_2
                    wgs_snv_file_pid2_original.vcf
    """
    # Read in the previously-created metadata dataframe.
    data = pd.read_csv(path_to_metadata_dataframe)
 
    # Create the initial folder (output_path_to_patient_folders) if not yet created. 
    if not os.path.exists(output_path_to_patient_folders):
        os.mkdir(output_path_to_patient_folders)

    # For each sample in the metadata dataframe, create a subfolder and copy the original SNV .vcf into the patient's subfolder.
    num_vcf_files_not_copied = 0
    for _, row in tqdm(data.iterrows(), total=data.shape[0], desc="Copying SNV Files"):
        pid = row["pid"]
        tumor_origin = row["tumor_origin"]
        path_to_vcf = row["path_to_wgs_file"]
        
        # If the .vcf file does not exist, skip.
        if not os.path.exists(path_to_vcf):
            num_vcf_files_not_copied += 1
            print(f"{path_to_vcf} does not exist")
            continue
        
        # Create folder for single patient if this does not exist.
        path_to_single_patient_folder = os.path.join(output_path_to_patient_folders, f"{pid}_{tumor_origin}")
        if not os.path.exists(path_to_single_patient_folder):
            os.mkdir(path_to_single_patient_folder)
            
        # If the patient's original .vcf file has not been copied 
        # to the patient subfolder, save .vcf file.
        name_of_new_vcf_file = os.path.join(
            path_to_single_patient_folder, 
            os.path.basename(path_to_vcf).replace(".vcf", "_original.vcf")
        )
        if not os.path.exists(name_of_new_vcf_file):
            shutil.copy(path_to_vcf, name_of_new_vcf_file)

    print(f".... Folder structure completed: {output_path_to_patient_folders} ")
    print(f".... Number of .vcf files not copied: {num_vcf_files_not_copied}")

def create_json_file(
    path_to_metadata_dataframe: str,
    path_to_patient_folders: str,
    output_path_to_results_json_file: str,
):
    """Create a .json file that will be used to track the results of the pipeline.

    Parameters
    ----------
    path_to_metadata_dataframe : str
        Path to the metadata dataframe.
    path_to_patient_folders : str
        Path to the folder containing all patient subfolders.
    output_path_to_results_json_file : str
        Output path to where you want to save the results .json file.
    """    
    
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
    
    metadata = pd.read_csv(path_to_metadata_dataframe)
    
    for folder in tqdm(glob.glob(os.path.join(path_to_patient_folders, "*"))):
        # Get the original .vcf filename.
        original_file = glob.glob(os.path.join(folder, f"*_original.vcf"))[0]
        
        # Get the PID from the .vcf file name.
        pid = _get_pid_from_structured_vcf_path(original_file, only_pid=True)
        assert pid in list(metadata["pid"]), f"{pid} does not exist in metadata dataframe."
        
        # Get the tumor_origin (primary tumor or metastasic tumor) and ge_data_available (True or False).
        tumor_origin = metadata[metadata["pid"] == pid].iloc[0]["tumor_origin"]

        if "path_to_rnaseq_file" in list(metadata.columns):
            ge_data_available = True if os.path.exists(metadata[metadata["pid"] == pid].iloc[0]["path_to_rnaseq_file"]) else False
        else:
            ge_data_available = metadata[metadata["pid"] == pid].iloc[0]["ge_data_available"]
            
        if (tumor_origin == "primary_tumor"):
            if (ge_data_available):
                results_dict["results"]["original"]["primary_tumor_wgs_and_rnaseq"].append(original_file)
            else:
                results_dict["results"]["original"]["primary_tumor_wgs"].append(original_file)
        else:
            if (ge_data_available):
                results_dict["results"]["original"]["metastasic_tumor_wgs_and_rnaseq"].append(original_file)
            else:
                results_dict["results"]["original"]["metastasic_tumor_wgs"].append(original_file)
                
    with open(output_path_to_results_json_file, "w") as f:
        json.dump(results_dict, f)
        
        