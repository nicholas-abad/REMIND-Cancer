import datetime
import json
import os
import pandas as pd
import subprocess
import time
from functools import wraps

def _wait_for_running_and_pending_lsf_cluster_jobs(
    maximum_number_of_jobs: int=0,
    sleep_timer: int=30
):
    """
    Wait for the number of running and pending jobs on an LSF cluster to drop below a specified maximum.
    Parameters
    ----------
    maximum_number_of_jobs : int, optional
        The maximum number of running and pending jobs to wait for. If this number is exceeded, the function will wait until the number drops below this value. Default is 0.
    sleep_timer : int, optional
        The number of seconds to wait between checks for the number of running and pending jobs. Default is 30.

    Returns
    -------
    None
    """
    assert type(maximum_number_of_jobs) == int
    assert type(sleep_timer) == int

    while True:
        # Get the number of running and pending jobs.
        num_running_and_pending_jobs = int(os.popen("bjobs -rp | wc -l").read())

        if num_running_and_pending_jobs <= maximum_number_of_jobs:
            break
        
        print(f"Waiting for all jobs to complete... Number of jobs: {num_running_and_pending_jobs}")
        time.sleep(sleep_timer)

    return

def _update_results_json_file(
    path_to_results_json: str,
    previous_tag_to_append: str,
    tag_to_append: str,
    dict_key: str,
    remaining_paths: list,
    filter_name: str,
):
    # Update the results .json file.
    with open(path_to_results_json, "r") as f:
        results = json.load(f)

    # Define the new new_tag within the dictionary.
    new_tag = previous_tag_to_append.replace(".vcf", tag_to_append)
    new_tag = previous_tag_to_append.replace(".vcf", tag_to_append)
    if previous_tag_to_append == "_original.vcf":
        previous_tag_to_append = "original"

    with open(path_to_results_json, "w") as f:
        # Check if the new_tag exists in the results dictionary. If not, create it.
        if new_tag not in results["results"]:
            results["results"][new_tag] = {}

        # Add the rows that still remain.
        results["results"][new_tag][dict_key] = remaining_paths

        # Overwrite the results_dict to include the new new_tag/value.
        json.dump(results, f)
        print(f"    {dict_key}:")
        print(
            f"        Total Number of Paths: {len(results['results'][previous_tag_to_append][dict_key])}"
        )
        print(f"        Number of Remaining Paths: {len(remaining_paths)}")


def remove_nested_parens(input_str: str) -> str:
    """Return a copy of 'input_str' with any parenthesized text removed.
    Nested parentheses are handled.
    
    Args:
    input_str (str): The string to remove nested parentheses from.
    
    Returns:
    str: A copy of 'input_str' with any parenthesized text removed.
    """
    result = []
    paren_level = 0
    for ch in input_str:
        if ch == "(":
            paren_level += 1
        elif ch == ")":
            paren_level -= 1
        else:
            if not paren_level:
                result.append(ch)
    return ''.join(result)

def _update_timing_tracker(
    path_to_vcf_file: str,
    name_of_pipeline_step: str,
    start_time: datetime,
    end_time: datetime,
    name_of_timing_tracker_file: str = "timing_tracker.json",
):
    # Ensure that name_of_timing_tracker_file ends with .json.
    if not name_of_timing_tracker_file.endswith(".json"):
        name_of_timing_tracker_file += ".json"
    
    # Get the path to the timing tracker.
    path_to_timing_tracker = os.path.join(
        os.path.dirname(path_to_vcf_file),
        name_of_timing_tracker_file
    )

    # If the timing tracker does not exist, create it as an empty dictionary.
    if not os.path.exists(path_to_timing_tracker):
        with open(path_to_timing_tracker, "w") as f:
            json.dump({}, f)
            
    # Load in the timing tracker.
    with open(path_to_timing_tracker, "r") as f:
        timing_tracker = json.load(f)
        
    # Create a key for the name of the pipeline step.
    if name_of_pipeline_step not in timing_tracker:
        timing_tracker[name_of_pipeline_step] = {
            "start_time": [],
            "end_time": [],
            "difference": []
        }
        
    # Calculate the difference and append these values to the dictionary.
    timing_tracker[name_of_pipeline_step]["start_time"].append(start_time.strftime("%d-%m-%Y_%H:%M"))
    timing_tracker[name_of_pipeline_step]["end_time"].append(end_time.strftime("%d-%m-%Y_%H:%M"))
    
    difference = end_time - start_time
    timing_tracker[name_of_pipeline_step]["difference"].append(difference.total_seconds())
    
    # Overwrite the time tracker.
    with open(path_to_timing_tracker, "w") as f:
        json.dump(timing_tracker, f)
    
            
def _get_pid_from_structured_vcf_path(
    path_to_vcf: str,
    only_pid: bool=False
):
    # Ex: H021-ABC_tumor (only_pid = False)
    # Ex: H021_ABC (only_pid = True)
    return path_to_vcf.split("/")[-2] if not only_pid else path_to_vcf.split("/")[-2].split("_")[0]

def execution_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Execution of '{func.__name__}' took {end_time - start_time:.4f} seconds.")
        return result
    return wrapper

def _get_index(header, column_name_desired, strict=False):
    header_split = header.split("\t")
    idx = -1
    for idx, col in enumerate(header_split):
        if not strict:
            if column_name_desired in col.rstrip():
                return idx

        # Consider two column names "TPM" and "z_score_TPM" and we're only looking for the "TPM" column.
        # If not strict, both would satisfy the above if-statement. This removes ambiguity.
        else:
            if column_name_desired == col.rstrip():
                return idx
    if idx == -1:
        return None


def _get_previous_tag(config: dict):
    with open(config["pipeline"]["path_to_results"], "r") as f:
        results = json.load(f)

    # Get the previous key (i.e. _after_promoter_filter.vcf)
    last_filter = list(results["results"].keys())[-1]

    # If this is the first filter, return an empty string.
    if last_filter == "original":
        last_filter = "_original.vcf"

    return last_filter

def _get_paths_from_previous_step(path_to_config: str):
    # Load in the configuration file.
    with open(path_to_config, "r") as f:
        config = json.load(f)

    # Get the SNV/Metastasis paths to consider.
    with open(config["pipeline"]["path_to_results"], "r") as f:
        # Load the dictionary that was initially created in Overview Patients (and continuously updated).
        paths_dict = json.load(f)
        # Get the previous step in the pipeline (i.e. original, after_cgc_filter, etc.)
        previous_step = list(paths_dict["results"].keys())[-1]

        # Get the four relevant paths.
        paths_to_tumor_snv_data = paths_dict["results"][previous_step]["tumor_snv"]
        paths_to_tumor_snv_and_ge_data = paths_dict["results"][previous_step][
            "tumor_snv_and_ge"
        ]
        paths_to_metastasis_snv_data = paths_dict["results"][previous_step][
            "metastasis_snv"
        ]
        paths_to_metastasis_snv_and_ge_data = paths_dict["results"][previous_step][
            "metastasis_snv_and_ge"
        ]

    return (
        paths_to_tumor_snv_data,
        paths_to_tumor_snv_and_ge_data,
        paths_to_metastasis_snv_data,
        paths_to_metastasis_snv_and_ge_data,
    )

def _create_cohort_dictionary(
    path_to_metadata: str,
):
    metadata = pd.read_csv(path_to_metadata)
    
    patient_cohort_dictionary = {}
    
    for idx, row in metadata.iterrows():
        pid = row["pid"]
        cohort = row["cohort"]
        if pid not in patient_cohort_dictionary:
            patient_cohort_dictionary[pid] = cohort
    return patient_cohort_dictionary

def _submit_bsub_job_and_get_job_id(
    bsub_pid: str="test",
    bsub_queue: str="short",
    bsub_allowed_time: str="0:10",
    bsub_memory_gb: str="1",
    command: str="ls",
):
    command = f'bsub -J {bsub_pid} -q {bsub_queue} -W {bsub_allowed_time} -M {bsub_memory_gb}GB -R rusage[mem={bsub_memory_gb}GB] "{command}"'
    job_id = int(subprocess.check_output(command, shell=True).decode("utf-8").strip().split()[1].replace("<", "").replace(">", ""))
    return job_id

def _get_bsub_job_statuses(job_ids):
    """
    Checks the status of a list of jobs and returns the summary status.

    Args:
        job_ids (list): A list of job IDs to check.

    Returns:
        dict: Returns a dictionary with the job IDs in each status category.
    """
    status_counts = {"DONE": [], "PEND": [], "EXIT": [], "UNKNOWN": []}

    for job_id in job_ids:
        try:
            job_status = subprocess.check_output(f"bjobs {job_id}", shell=True).decode("utf-8").strip().split("\n")[-1]
        except subprocess.CalledProcessError:
            # If the `bjobs` command fails, consider the status as "UNKNOWN".
            status_counts["UNKNOWN"].append(job_id)
        else:
            if "DONE" in job_status:
                status_counts["DONE"].append(job_id)
            elif "PEND" in job_status:
                status_counts["PEND"].append(job_id)
            elif "EXIT" in job_status:
                status_counts["EXIT"].append(job_id)
            else:
                status_counts["UNKNOWN"].append(job_id)

    return status_counts

def _wait_for_bsub_jobs_completion(job_ids, wait_time):
    """
    Waits until all job IDs are in the "DONE" status. If any job IDs are in the "EXIT" status, raises an error.

    Args:
        job_ids (list): A list of job IDs to wait for.

    Raises:
        RuntimeError: If one or more job IDs are in the "EXIT" status.
    """
    while True:
        status_counts = _get_bsub_job_statuses(job_ids)

        if len(status_counts["EXIT"]) > 0:
            # If any job IDs are in the "EXIT" status, raise an error.
            raise RuntimeError("One or more job IDs are in the 'EXIT' status.")

        if len(status_counts["DONE"]) == len(job_ids):
            # If all job IDs are in the "DONE" status, break out of the loop.
            print(f"All {len(job_ids)} job_ids are complete!")
            break

        # Wait for wait_time seconds before checking the status again.
        print(f"DONE: {len(status_counts['DONE'])}/{len(job_ids)}")
        time.sleep(wait_time)

def _get_number_of_lines(path, include_header_in_count=False):
    result = subprocess.run(['wc', '-l', path], capture_output=True, text=True)
    num_lines = int(result.stdout.strip().split()[0])
    if not include_header_in_count:
        num_lines -= 1
    return num_lines

def parse_date(sequencing_date: str):
    """Parse the date into the proper format for easy comparison.

    Four possible formats that the sequencing date comes in:
        - "2017-05-15_13h21"
        - "2018-05-15_13h21_+0200"
        - "2018-05-15_13h21_Europe_Berlin"
        - ""

    This takes the first 17 characters (i.e. up to 2017-05-15_13h21) and parses it into a datetime object.

    Parameters
    ----------
    sequencing_date : str
        Sequencing date given as a string.

    Returns
    -------
    datetime object
        A datetime object that corresponds to the given sequencing date.
    """
    # If empty, give it a random early date.
    if sequencing_date == "":
        return datetime.datetime.strptime("2001-05-15_13h21", "%Y-%m-%d_%Hh%M")

    else:
        sequencing_date = sequencing_date[:16]
        dt = datetime.datetime.strptime(sequencing_date, "%Y-%m-%d_%Hh%M")
        return dt


def _get_most_recent_vcf_file(line1, line2):
    date_of_line1 = parse_date("_".join(line1.split("/")[-2].split("_")[4:]))
    date_of_line2 = parse_date("_".join(line2.split("/")[-2].split("_")[4:]))

    # Date1 occurs later than date2.
    if date_of_line1 > date_of_line2:
        return line1
    else:
        return line2
    
def _search_for_gene_in_pid(
    pid_to_search_for: str,
    gene_to_search_for: str,
    results: dict
):

    information = {}

    for filter_step in results["results"]:
        pid_in_filter_step = False
        gene_in_filter_step = False
        path_of_pid = None

        if pid_to_search_for not in information:
            information[pid_to_search_for] = {}

        for category in results["results"][filter_step]:
            for path in results["results"][filter_step][category]:
                if pid_to_search_for in path:
                    pid_in_filter_step = True
                    filter_step_df = pd.read_csv(path, usecols=["#CHROM", "GENE", "REF", "ALT"], delimiter="\t")
                    gene_in_filter_step = True if gene_to_search_for in list(filter_step_df["GENE"]) else False
                    path_of_pid = path

        information[pid_to_search_for][filter_step] = {
            "pid_in_filter_step_bool": pid_in_filter_step,
            "gene_in_filter_step": gene_in_filter_step,
            "path": path_of_pid
        }

    print(json.dumps(information, indent=4))
    
def _get_statistics_per_filter_for_pid(
    pid_to_search_for: str,
    results: dict,
):
    """
    For a pid, check whether they pass each filter. 
    Also, check for the number of lines in each filter.
    """
    information = {}

    for filter_step in results["results"]:
        pid_in_filter_step = False
        path_of_pid = None
        number_of_lines = None

        if pid_to_search_for not in information:
            information[pid_to_search_for] = {}

        for category in results["results"][filter_step]:
            for path in results["results"][filter_step][category]:
                if pid_to_search_for in path:
                    pid_in_filter_step = True
                    path_of_pid = path
                    number_of_lines = _get_number_of_lines(path)
        information[pid_to_search_for][filter_step] = {
            "pid_in_filter_step_bool": pid_in_filter_step,
            "path_of_pid": path_of_pid,
            "number_of_lines": number_of_lines
        }
    print(json.dumps(information, indent=4))