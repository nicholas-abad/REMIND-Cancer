import json
import os
from tqdm import tqdm
import subprocess

from optparse import OptionParser

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)
    
from src.pipeline.general_helper_functions import (
    _get_pid_from_structured_vcf_path,
    _update_results_json_file, _wait_for_running_and_pending_lsf_cluster_jobs,
    _submit_bsub_job_and_get_job_id
    )

def main(path_to_config: str):
    # Load in the configuration file.
    with open(path_to_config, "r") as f:
        config = json.load(f)
        
    # Get the path to the single file python script.
    path_to_single_file_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "_run_ge_on_single_file.py",
    )
    
    assert os.path.exists(path_to_single_file_script)
    
    run_on_cluster = config["pipeline"]["run_on_cluster"]
    
    with open(config["pipeline"]["path_to_results"], "r") as f:
        results = json.load(f)

    # Get the previous step in the pipeline (i.e. original, after_cgc_filter, etc.)
    previous_step = list(results["results"].keys())[-1]
    suffix_to_append_to_vcf = config["ge_filter"]["suffix_to_append_to_vcf"]


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

                job_id = _submit_bsub_job_and_get_job_id(
                    bsub_pid = pid,
                    bsub_queue = "medium",
                    bsub_allowed_time = "1:00",
                    bsub_memory_gb = "5",
                    command = command,
                )
                # job_ids_to_wait_for.append(job_id)

            else:
                # Use Popen to execute the command and stream the output in real-time
                with subprocess.Popen(command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
                    for line in proc.stdout:
                        print(line, end='')  # Print the output in real-time
                    for err in proc.stderr:
                        print(err, end='')   # Print any error messages in real-time
        
        # _wait_for_bsub_jobs_completion(job_ids_to_wait_for, wait_time=5)
    
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
            previous_tag_to_append=previous_step,
            tag_to_append=suffix_to_append_to_vcf,
            dict_key=key,
            remaining_paths=remaining_paths,
            filter_name="ge_filter"
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
    
    