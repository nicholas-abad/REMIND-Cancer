import pandas as pd
from tqdm import tqdm
from optparse import OptionParser
import json
import numpy as np
import os
from datetime import datetime

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)
    
from src.pipeline.general_helper_functions import (execution_time, _get_index,
                                                   _get_pid_from_structured_vcf_path,
                                                   _create_cohort_dictionary, _update_timing_tracker)


def _create_recurrence_dict(
    paths: list,
    expression_score_to_use: str,
    cohort_dictionary: dict,
):
    """
    Iterate through all original paths and create a dictionary based off of their mutation location.
    The position is the key and the value of that key is a list that contains all other mutations within the
    entire dataset that also have that basepair location. Each entry in that list is a tuple that contains
    (1) the filepath of the mutation, (2) the patient id, (3) the reference nucleotide, and (4) the
    alternate nucleotide.

    """

    # Gather all of the base positions.
    recurrence_dict = {}
    """
    {
        base_position: [
            (filepath, pid, reference_nucleotide, alt_nucleotide),
            (filepath, pid, reference_nucleotide, alt_nucleotide)
            ]
    }
    """
    files_that_do_not_exist = []
    for file in tqdm(paths):
        with open(file, "r") as f:
            header = f.readline()

            # Get indices of relevant columns
            base_position_idx = _get_index(header, "POS", True)
            reference_idx = _get_index(header, "REF", True)
            alternate_idx = _get_index(header, "ALT", True)
            gene_idx = _get_index(header, "GENE", True)
            chromosome_idx = _get_index(header, "#CHROM", True)

            confidence_idx = _get_index(header, "CONFIDENCE", True)
            raw_score_idx = _get_index(header, expression_score_to_use, True)
            zscore_idx = _get_index(
                header, f"{expression_score_to_use}_Z_score", True)

            purity_median_idx = _get_index(header, "purity", True)
            allele_frequency_idx = _get_index(header, "allele_frequency", True)

            for line in f:
                line_split = line.replace("\n", "").split("\t")

                # Get the current bp, reference, alternative and chromosome values.
                current_base_position = line_split[base_position_idx]
                current_ref = line_split[reference_idx]
                current_alt = line_split[alternate_idx]
                current_gene = line_split[gene_idx]
                current_chromosome = str(line_split[chromosome_idx])

                if raw_score_idx != None:
                    current_raw_score = line_split[raw_score_idx] if line_split[raw_score_idx] != "" else 0
                    current_zscore = line_split[zscore_idx] if line_split[raw_score_idx] != "" else 0
                    current_log_score = np.log(float(
                        current_raw_score)) if current_raw_score != "not_available" else "not_available"
                else:
                    current_raw_score = "not_available"
                    current_zscore = "not_available"
                    current_log_score = "not_available"

                current_confidence = line_split[confidence_idx] if confidence_idx != None else "not_available"
                current_purity = line_split[purity_median_idx] if purity_median_idx != None else "not_available"
                current_allele_frequency = line_split[
                    allele_frequency_idx] if allele_frequency_idx != None else "not_available"

                if current_chromosome not in recurrence_dict:
                    recurrence_dict[current_chromosome] = {}

                if current_base_position not in recurrence_dict[current_chromosome]:
                    recurrence_dict[current_chromosome][current_base_position] = [
                    ]

                patient_id = _get_pid_from_structured_vcf_path(
                    file, only_pid=True)
                patient_cohort = cohort_dictionary[patient_id] if patient_id in cohort_dictionary else "not_available"

                recurrence_dict[current_chromosome][current_base_position].append(
                    (file, patient_id, patient_cohort, current_base_position, current_ref, current_alt, current_gene, current_chromosome,
                     current_raw_score, current_zscore, current_log_score, current_confidence, current_purity, current_allele_frequency)
                )

    return recurrence_dict, files_that_do_not_exist


@execution_time
def _add_recurrence(
    path_to_vcf_file: str,
    compute_recurrence_dictionary_on_current_dataset: bool,
    additional_recurrence_datasets_to_add: list,
    name_of_column_to_add: str,
    dataset: str,
    ending_of_filter_step: str,
    rnaseq_measurement: str,
    results: dict,
    suffix_to_append_to_vcf: str,
    path_to_metadata: str,
):
    # Load in the dataframe from the .vcf file.
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")
    pid = _get_pid_from_structured_vcf_path(path_to_vcf_file, True)

    # Get the cohort of the current patient.
    patient_cohort_dictionary = _create_cohort_dictionary(path_to_metadata)
    patient_cohort = patient_cohort_dictionary[pid]

    # Define a dictionary such that each key is the dataset name and the value is that datasets recurrence dictionary.
    recurrence_dictionaries_to_consider = {}

    # Iterate through all additional_recurrence_datasets_to_add and add this to recurrence_dictionaries_to_consider.
    for dataset_name, path_to_recurrence_dictionary in additional_recurrence_datasets_to_add:
        if os.path.exists(path_to_recurrence_dictionary):
            print(
                f"Loading {dataset_name} dataset for recurrence ({path_to_recurrence_dictionary})")
            assert dataset_name not in recurrence_dictionaries_to_consider, f"Name of recurrence dictionaries need to be unique!"
            with open(path_to_recurrence_dictionary, "r") as f:
                recurrence_dictionaries_to_consider[dataset_name] = json.load(
                    f)

    # If compute_recurrence_dictionary_on_current_dataset, add this dataset as well as the recurrence dictionary to the list as well.
    if compute_recurrence_dictionary_on_current_dataset:
        # To compute the recurrence dictionary, get the ending of the paths to consider.
        # By default, only files after the promoter filter are considered.
        results_key_after_promoter_filter = [key for key in list(
            results["results"].keys()) if key.endswith(ending_of_filter_step)][0]

        paths_to_consider = []
        for category in results["results"][results_key_after_promoter_filter]:
            paths_to_consider += results["results"][results_key_after_promoter_filter][category]

        print(
            f"Calculating recurrence dictionary for {len(paths_to_consider)} paths.")
        print(
            f"  Paths end with the suffix: '{results_key_after_promoter_filter}'")

        recurrence_dict, _ = _create_recurrence_dict(
            paths_to_consider, rnaseq_measurement, patient_cohort_dictionary)
        recurrence_dictionaries_to_consider[f"current_{dataset}_dataset"] = recurrence_dict

    # Iterate through the recurrence dictionaries and add the recurrence to the dataframe.
    for name_of_recurrence_dictionary in recurrence_dictionaries_to_consider:
        recurrence_dictionary_to_consider = recurrence_dictionaries_to_consider[
            name_of_recurrence_dictionary]

        for idx, row in data.iterrows():
            chromosome = str(row["#CHROM"])
            position = str(row["POS"])
            ref = row["REF"]
            alt = row["ALT"]
            gene = row["GENE"]

            recurrence_column_value = ""
            recurrence_count = 0

            if chromosome in recurrence_dictionary_to_consider:
                if position in recurrence_dictionary_to_consider[chromosome]:
                    for possible_recurrence_entry in recurrence_dictionary_to_consider[chromosome][position]:
                        rec_path, rec_pid, rec_cohort, rec_bp, rec_ref, rec_alt, rec_gene, rec_chromosome, rec_exp, rec_zscore, rec_log, \
                            rec_confidence, rec_purity, rec_af = possible_recurrence_entry
                        if (pid != rec_pid) and (alt == rec_alt) and (gene == rec_gene):
                            recurrence_column_value += ",".join(
                                [str(i) for i in possible_recurrence_entry]) + ";"
                            recurrence_count += 1
            data.loc[idx, f"{name_of_recurrence_dictionary}_{name_of_column_to_add}"] = recurrence_column_value[:-
                                                                                                                1] if len(recurrence_column_value) != 0 else recurrence_column_value
            data.loc[idx, f"{name_of_recurrence_dictionary}_num_recurrent_mutations"] = recurrence_count

    # Write the .csv.
    data.to_csv(
        path_to_vcf_file.replace(".vcf", suffix_to_append_to_vcf),
        sep="\t", index=False
    )


@execution_time
def main(
    path_to_vcf_file: str,
    config: dict
):
    compute_recurrence_dictionary_on_current_dataset = config[
        "recurrence"]["compute_recurrence_with_current_dataset"]
    additional_recurrence_datasets_to_add = config["recurrence"]["additional_recurrence_datasets_to_add"]
    name_of_column_to_add = config["recurrence"]["name_of_column_to_add"]
    suffix_to_append_to_vcf = config["recurrence"]["suffix_to_append_to_vcf"]

    dataset = config['pipeline']['dataset']
    ending_of_filter_step = config["promoter"]["suffix_to_append_to_vcf"]
    rnaseq_measurement = config["rnaseq_data"]["rnaseq_measurement"]
    path_to_results = config["pipeline"]["path_to_results"]
    path_to_metadata = config["pipeline"]["path_to_metadata"]

    with open(path_to_results, "r") as f:
        results = json.load(f)

    _add_recurrence(
        path_to_vcf_file=path_to_vcf_file,
        compute_recurrence_dictionary_on_current_dataset=compute_recurrence_dictionary_on_current_dataset,
        additional_recurrence_datasets_to_add=additional_recurrence_datasets_to_add,
        name_of_column_to_add=name_of_column_to_add,
        dataset=dataset,
        ending_of_filter_step=ending_of_filter_step,
        rnaseq_measurement=rnaseq_measurement,
        results=results,
        suffix_to_append_to_vcf=suffix_to_append_to_vcf,
        path_to_metadata=path_to_metadata,
    )


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--path-to-vcf-file",
        action="store",
        type="str",
        dest="path_to_vcf_file",
    )
    parser.add_option(
        "--config",
        action="store",
        type="str",
        dest="config",
    )

    (options, args) = parser.parse_args()

    path_to_vcf_file = options.path_to_vcf_file
    path_to_config = options.config

    with open(path_to_config, "r") as f:
        config = json.load(f)

    start_time = datetime.now()

    main(
        path_to_vcf_file=path_to_vcf_file,
        config=config
    )
    
    end_time = datetime.now()
    
    _update_timing_tracker(
        path_to_vcf_file = path_to_vcf_file,
        name_of_pipeline_step = "recurrence",
        start_time = start_time,
        end_time = end_time,
        name_of_timing_tracker_file = "timing_tracker.json",
    )

    print("~~ Finished with adding recurrence ~~")
