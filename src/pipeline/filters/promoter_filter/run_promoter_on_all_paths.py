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
    _submit_bsub_job_and_get_job_id
    )


def main(path_to_config: str):
    # Load in the configuration file.
    with open(path_to_config, "r") as f:
        config = json.load(f)
        
    # Get the path to the single file python script.
    path_to_single_file_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "_run_promoter_on_single_file.py",
    )
    
    assert os.path.exists(path_to_single_file_script)
    
    run_on_cluster = config["pipeline"]["run_on_cluster"]
    
    with open(config["pipeline"]["path_to_results"], "r") as f:
        results = json.load(f)

    # Get the previous step in the pipeline (i.e. original, after_cgc_filter, etc.)
    previous_step = list(results["results"].keys())[-1]
    suffix_to_append_to_vcf = config["promoter"]["suffix_to_append_to_vcf"]

    # Get the four relevant lists of paths.
    paths_to_tumor_snv_data = results["results"][previous_step]["tumor_snv"]
    paths_to_tumor_snv_and_ge_data = results["results"][previous_step]["tumor_snv_and_ge"]
    paths_to_metastasis_snv_data = results["results"][previous_step]["metastasis_snv"]
    paths_to_metastasis_snv_and_ge_data = results["results"][previous_step]["metastasis_snv_and_ge"]
    
    paths_and_keys = (
        (paths_to_tumor_snv_data, "tumor_snv"),
        (paths_to_tumor_snv_and_ge_data, "tumor_snv_and_ge"),
        (paths_to_metastasis_snv_data, "metastasis_snv"),
        (paths_to_metastasis_snv_and_ge_data, "metastasis_snv_and_ge"),
    )
        
    for paths, key in paths_and_keys:
        # job_ids_to_wait_for = []
        for path in tqdm(paths, total=len(paths), desc=key):
            suffix_to_append_to_vcf = config["promoter"]["suffix_to_append_to_vcf"]

            command = f"python {path_to_single_file_script} --path-to-vcf-file {path} --config {path_to_config}"

            if run_on_cluster:
                _wait_for_running_and_pending_lsf_cluster_jobs(
                    maximum_number_of_jobs = 300,
                    sleep_timer = 10
                )

                pid = _get_pid_from_structured_vcf_path(path, True)

                pids_with_resource_limitations = [
                    # NCT-MASTER (retrospective)
                    "H021-QRG2SP",
                    "H021-JRJTPJ",
                    "H021-V48HPK",
                    "H021-GZV3MZ",
                    "H021-4992WJ",
                    'H021-LDDQVD',
                    'H021-JRJTPJ',
                    'H021-DTEPSQ',
                    'H021-WC2X6E',
                    'H021-5KYMBG',
                    'H021-FXCFXC',
                    'H021-64F9KW',
                    'H021-AR6ZSW',
                    'H021-5BJNYD',
                    'H021-GKBJ69',
                    'H021-MY6WP7',
                    'H021-HBTTLZ',
                    'H021-YUKAHB',
                    'H021-6ERRP7',
                    'H021-T5UEDU',
                    'H021-FS6HNA',
                    'H021-737YQZ',
                    'H021-MKVVUC',
                    'H021-V48HPK',
                    'H021-FRJ1TB',
                    'H021-DMJ8',
                    'H021-9XB896',
                    'H021-Y78RB8',
                    'H021-NVGQD4',
                    'H021-PXVWR6',
                    'H021-LJMFK2',
                    'H021-UA82VZ',
                    'H021-9GCB23',
                    'H021-HVX6T8',
                    'H021-C41MJC',
                    'H021-D6925N',
                    'H021-RF6KBF',
                    'H021-6181DL',
                    'H021-4992WJ',
                    'H021-ZSW8VM',
                    'H021-QRG2SP',
                    'H021-DKK8FG',
                    'H021-9QZTUK',
                    'H021-7HFAWC',
                    'H021-SGXVMY',
                    'H021-SDXUG8',
                    'H021-NSQFWT',
                    'H021-U67T4B',
                    'H021-GZV3MZ',
                    'H021-YHPKYT',
                    'H021-2MJBY3',
                    'H021-H4QNNU',
                    'H021-AJAAXM',
                    'H021-6F52XV',
                    'H021-89D53Z',
                    'H021-DTHS9M',
                    'H021-XSCYZ5',
                    'H021-RKGWNV',
                    'H021-Y8G3M6',
                    'H021-MTNLKM',
                    'H021-9QCLQ4',
                    'H021-VB46RR',
                    'H021-25LC94',
                    'H021-EBYSD9',
                    'H021-XJTUSV',
                    # PCAWG
                    "04aa6b77-8074-480c-872e-a1a47afa5314", # Done
                    "0980e7fd-051d-45e9-9ca6-2baf073da4e8", # Done
                    "c8e961b4-e324-40a2-89f6-736ec3845bc9", # Done
                    "93ff786e-0165-4b02-8d27-806d422e93fc", # Done
                    "3869ff3f-21b9-4817-8ff4-83c6fc75ab11", # Done
                    "75ba6722-1148-4a52-a9ed-68d890238205", # Done
                    "6ca5c1bb-275b-4d05-948a-3c6c7d03fab9", # Done
                    "ed32c725-08ae-48eb-8fa2-719b9aeb7550", # Length mismatch > Fixed
                    "f92a34fa-014e-4b41-a6d0-3b46b8c8a3ee", # Length mismatch > Fixed
                    "3f98d326-5676-4257-9af8-0a5f5d3c2527", # Length mismatch > Fixed
                    "e7d74d34-3255-4c20-90fd-b105e6e229c8", # Length mismatch > Fixed
                    "14c5b81d-da49-4db1-9834-77711c2b1d38", # Done
                    "2df02f2b-9f1c-4249-b3b4-b03079cd97d9", # Done
                    "5c3def3a-b515-41f6-8157-681b963534e7", # Done
                    "98e8f23c-5970-4fce-9551-4b11a772fe1b", # Done
                    "deb9fbb6-656b-41ce-8299-554efc2379bd", # Done
                    "60413de1-6cd2-4f74-8180-3bdd394d6d16", # Done
                    "63762458-902a-4329-a823-703b54cb5f9d", # Done
                    "2790b964-63e3-49aa-bf8c-9a00d3448c25", # Done
                    "51893d3f-e7f3-43f9-9fd0-c0f25ae96804", # Done
                    "ca8fa9f5-3190-440d-9879-22e33d05ca6c", # Done
                    "418e916b-7a4e-4fab-8616-15dcec4d79f8", # Length mismatch > Fixed
                ]
                
                
                if pid in pids_with_resource_limitations:
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
                # job_ids_to_wait_for.append(job_id)

            else:
                os.popen(command)
        
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
            filter_name="promoter"
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
    
    