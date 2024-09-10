import json
import os
from tqdm import tqdm

from optparse import OptionParser

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import (
    _get_pid_from_structured_vcf_path,
    _update_results_json_file, _wait_for_running_and_pending_lsf_cluster_jobs,
    _submit_bsub_job_and_get_job_id, _get_number_of_lines
    )


def main(path_to_config: str):
    # Load in the configuration file.
    with open(path_to_config, "r") as f:
        config = json.load(f)
        
    # Get the path to the single file python script.
    path_to_single_file_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "_run_preprocessing_on_single_path.py",
    )
    assert os.path.exists(path_to_single_file_script)
    
    run_on_cluster = config["pipeline"]["run_on_cluster"]
    
    with open(config["pipeline"]["path_to_results"], "r") as f:
        results = json.load(f)

    # Get the previous step in the pipeline (i.e. original, after_cgc_filter, etc.)
    # previous_step = list(results["results"].keys())[-1]
    previous_step = "original"
    print("Previous Step: ", previous_step)
    #previous_suffix_to_append_to_vcf = _get_previous_tag(config)
    previous_suffix_to_append_to_vcf = "_original.vcf"
    suffix_to_append_to_vcf = config["preprocessing_details"]["suffix_to_append_to_vcf"]


    # Get the four relevant lists of paths.
    paths_to_tumor_snv_data = results["results"][previous_step]["primary_tumor_wgs"]
    paths_to_tumor_snv_and_ge_data = results["results"][previous_step]["primary_tumor_wgs_and_rnaseq"]
    paths_to_metastasis_snv_data = results["results"][previous_step]["metastasic_tumor_wgs"]
    paths_to_metastasis_snv_and_ge_data = results["results"][previous_step]["metastasic_tumor_wgs_and_rnaseq"]
    
    paths_and_keys = (
        (paths_to_tumor_snv_data, "primary_tumor_wgs"),
        (paths_to_tumor_snv_and_ge_data, "primary_tumor_wgs_and_rnaseq"),
        (paths_to_metastasis_snv_data, "metastasic_tumor_wgs"),
        (paths_to_metastasis_snv_and_ge_data, "metastasic_tumor_wgs_and_rnaseq"),
    )
        
    for paths, key in paths_and_keys:
        # job_ids_to_wait_for = []
        for path in tqdm(paths, total=len(paths), desc=key):

            command = f"python {path_to_single_file_script} --path-to-vcf-file {path} --config {path_to_config}"

            if run_on_cluster:
                _wait_for_running_and_pending_lsf_cluster_jobs(
                    maximum_number_of_jobs = 300,
                    sleep_timer = 10
                )

                pid = _get_pid_from_structured_vcf_path(path, True)
                
                # Check how many lines in the file in order to send to proper queue in cluster.
                num_lines = _get_number_of_lines(path)
                if num_lines >= 100000:                    
                    job_id = _submit_bsub_job_and_get_job_id(
                        bsub_pid = pid,
                        bsub_queue = "long",
                        bsub_allowed_time = "10:00",
                        bsub_memory_gb = "50",
                        command = command,
                    )
                else:
                    job_id = _submit_bsub_job_and_get_job_id(
                        bsub_pid = pid,
                        bsub_queue = "medium",
                        bsub_allowed_time = "1:00",
                        bsub_memory_gb = "5",
                        command = command,
                    )

            else:
                os.popen(command)
            
        if run_on_cluster:
            _wait_for_running_and_pending_lsf_cluster_jobs(
                maximum_number_of_jobs = 0,
                sleep_timer = 10
            )
        
        remaining_paths = []
        
        for original_path in paths:
            new_path = original_path.replace(".vcf", suffix_to_append_to_vcf)
            if os.path.exists(new_path):
                remaining_paths.append(new_path)

        _update_results_json_file(
            path_to_results_json=config["pipeline"]["path_to_results"],
            previous_tag_to_append=previous_suffix_to_append_to_vcf,
            tag_to_append=suffix_to_append_to_vcf,
            dict_key=key,
            remaining_paths=remaining_paths,
            filter_name="preprocessing"
        )
    
    

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--config",
        action="store",
        type="str",
        dest="config",
    )
        
    (options, args) = parser.parse_args()
    
    path_to_config = options.config
    
    main(path_to_config)         
    
    