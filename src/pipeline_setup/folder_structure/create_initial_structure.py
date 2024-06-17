from optparse import OptionParser
from _main_functions import create_new_structure, create_json_file
import json

def main(path_to_config: str):
    with open(path_to_config, "r") as f:
        config = json.load(f)

    # Define the path to where to save the pipeline results .json file.
    # NOTE: This .json file keeps track of results and .vcf filenames.
    output_path_to_results_json_file = config["pipeline"]["path_to_results"]

    # Path to the previously-created metadata dataframe.
    path_to_metadata_dataframe = config["pipeline"]["path_to_metadata"]

    # Name of the dataset.
    dataset = config["pipeline"]["dataset"]

    # Define the path to where to copy the original data files (located within the metadataframe)
    # to the folder structure that's applicable with the REMIND-Cancer Pipeline.
    output_path_to_patient_folders = config["pipeline_folder_setup"][dataset]["output_path_to_patient_folders"]
    
    # Define an extra suffix to add to indicate the "original" SNV files.
    # NOTE: The "original" SNV files refers to the copied version that will be within 'output_path_to_patient_folders' 
    ending_of_original_files = config["pipeline_folder_setup"][dataset]["ending_suffix_to_add_to_original_files"]

    # Create the folder structure.
    create_new_structure(
        path_to_metadata_dataframe = path_to_metadata_dataframe,
        path_to_patient_folders = output_path_to_patient_folders,
        ending_of_original_files = ending_of_original_files    
    )

    # Create the initial .json file that will be used throughout the pipeline.
    create_json_file(
        path_to_metadata_dataframe = path_to_metadata_dataframe,
        path_to_patient_folders = output_path_to_patient_folders,
        path_to_results_json_file = output_path_to_results_json_file,
        ending_of_original_files = ending_of_original_files,
    )
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--path-to-config",
        action="store",
        type="str",
        dest="config",
    )
        
    (options, args) = parser.parse_args()
    
    path_to_config = options.config
    
    main(path_to_config)    