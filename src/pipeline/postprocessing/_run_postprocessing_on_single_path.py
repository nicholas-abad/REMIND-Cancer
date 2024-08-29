import json
import os
import pandas as pd
from tqdm import tqdm

import numpy as np
from scipy.stats import zscore
from optparse import OptionParser

from datetime import datetime

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import (
    _get_number_of_lines, 
    _get_pid_from_structured_vcf_path, 
    _create_cohort_dictionary,
    _update_timing_tracker
)


def _add_tf_logo_sequence_to_single_patient_path(
    dataframe: pd.DataFrame,
    path_to_downloaded_jaspar_json: str,

) -> pd.DataFrame:

    # (1) Generate a dictionary for TF names to sequence logos.
    # Ex: tf_name_to_matrix_id_and_sequence_logos[tf_name]["sequence_logo"] = https://jaspar2022.genereg.net/static/logos/...

    tf_name_to_matrix_id_and_sequence_logos = {}
    assert os.path.exists(
        path_to_downloaded_jaspar_json), f"You first need to download the jaspar file: https://jaspar2022.genereg.net/api/v1/matrix/?format=json&page_size=1000&collection=CORE&tax_group=Vertebrates&tax_id=9606"

    with open(path_to_downloaded_jaspar_json, "r") as f:
        downloaded_jaspar_response = json.load(f)
        for tf_entry in tqdm(downloaded_jaspar_response["results"]):
            tf_name = tf_entry["name"]
            matrix_id = tf_entry["matrix_id"]
            sequence_logo = f"https://jaspar2022.genereg.net/static/logos/all/{matrix_id}.png"
            rc_sequence_logo = f"https://jaspar2022.genereg.net/static/logos/all/{matrix_id}.rc.png"
            version = int(tf_entry["version"])

            # Get the newer version.
            if tf_name in tf_name_to_matrix_id_and_sequence_logos:
                prev_version = tf_name_to_matrix_id_and_sequence_logos[tf_name]["version"]
                if version < prev_version:
                    continue

            tf_name_to_matrix_id_and_sequence_logos[tf_name] = {
                "matrix_id": matrix_id,
                "sequence_logo": sequence_logo,
                "rc_sequence_logo": rc_sequence_logo,
                "version": version
            }

    # (2) Add sequence logo to the patient.
    name_of_jaspar_column = config["preprocessing_details"]["transcription_factor_prediction"][
        "name_of_tfbs_prediction_column_to_add"] + "(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"
    new_name_of_jaspar_column = config["preprocessing_details"]["transcription_factor_prediction"][
        "name_of_tfbs_prediction_column_to_add"] + "(tf_name,binding_affinity,seq1,seq2,raw,zscore,log,tf_sequence_logo)"

    assert name_of_jaspar_column in list(dataframe.columns)

    for idx, row in dataframe.iterrows():
        jaspar_entries = row[name_of_jaspar_column]
        if jaspar_entries == ".":
            dataframe.loc[idx, new_name_of_jaspar_column] = "."
        else:
            new_jaspar_entries = []
            for entry in jaspar_entries.split(";"):
                entry_split = entry.split(",")
                tf_name = entry_split[0]
                if tf_name in tf_name_to_matrix_id_and_sequence_logos:
                    entry_split.append(
                        tf_name_to_matrix_id_and_sequence_logos[tf_name]["sequence_logo"])
                else:
                    entry_split.append("not_available")
                combined_entry_split = ",".join(entry_split)
                new_jaspar_entries.append(combined_entry_split)
            dataframe.loc[idx, new_name_of_jaspar_column] = ";".join(
                new_jaspar_entries)

    return dataframe


def _add_num_original_promoter_and_final_mutations_to_single_patient_path(
    dataframe: pd.DataFrame,
    pid: str,
    results: dict,
    ending_after_preprocessing: str,
    ending_after_promoter: str,
):
    # Get the original (after preprocessing) and promoter results key.
    list_of_results_keys = [key for key in results["results"].keys()]

    after_preprocessing_results_key = [
        key for key in list_of_results_keys if key.endswith(ending_after_preprocessing)][0]
    after_promoter_results_key = [
        key for key in list_of_results_keys if key.endswith(ending_after_promoter)][0]

    # Get the number of mutations after preprocessing.
    for category in results["results"][after_preprocessing_results_key]:
        for file in results["results"][after_preprocessing_results_key][category]:
            if pid in file:
                num_mutations_after_preprocessing = _get_number_of_lines(file)
                break

    # Get the number of mutations after the promoter filter.
    for category in results["results"][after_promoter_results_key]:
        for file in results["results"][after_promoter_results_key][category]:
            if pid in file:
                num_mutations_after_promoter_filter = _get_number_of_lines(
                    file)
                break

    # Get the number of mutations in this current file.
    num_mutations_after_pipeline = dataframe.shape[0]

    dataframe["num_original_mutations"] = num_mutations_after_preprocessing
    dataframe["num_promoter_mutations"] = num_mutations_after_promoter_filter
    dataframe["num_final_mutations"] = num_mutations_after_pipeline

    return dataframe


def _add_ncbi_information_of_genes_and_tfs(
    dataframe: pd.DataFrame,
    path_to_gene_name_and_description_file: str,
    name_of_tfbs_prediction_column: str,
):
    with open(path_to_gene_name_and_description_file, "r") as f:
        gene_name_and_description_file = json.load(f)

    for idx, row in dataframe.iterrows():
        gene = row["GENE"]
        tf_column = row[name_of_tfbs_prediction_column]

        list_of_genes_and_tfs = [gene]

        if tf_column != ".":
            for entry in tf_column.split(";"):
                tf_name = entry.split(",")[0]
                list_of_genes_and_tfs.append(tf_name)

        gene_and_tf_description_single_row = {}

        for gene_or_tf in list_of_genes_and_tfs:
            if gene_or_tf in gene_name_and_description_file:
                description = gene_name_and_description_file[gene_or_tf]
            else:
                description = "unknown"
            gene_and_tf_description_single_row[gene_or_tf] = description
        dataframe.loc[idx, "ncbi_gene_and_tf_summaries"] = str(
            gene_and_tf_description_single_row)

    return dataframe


def _get_pid_from_structured_vcf_path(
    path_to_vcf: str,
    only_pid: bool = False
):
    # Ex: H021-ABC_tumor (only_pid = False)
    # Ex: H021_ABC (only_pid = True)
    return path_to_vcf.split("/")[-2] if not only_pid else path_to_vcf.split("/")[-2].split("_")[0]


def _get_relevant_cohorts_for_single_row(row: pd.Series, cohort_of_patient: str, name_of_recurrence_column: str):
    relevant_cohorts = [cohort_of_patient]
    recurrence_entries = row[name_of_recurrence_column]
    if pd.isna(recurrence_entries) or recurrence_entries == ".":
        return relevant_cohorts
    for entry in recurrence_entries.split(";"):
        cohort = entry.split(",")[2]
        relevant_cohorts.append(cohort)
    relevant_cohorts = list(set(relevant_cohorts))
    return relevant_cohorts


def _get_relevant_genes_and_tfs_for_single_row(row: pd.Series, current_gene: str, name_of_jaspar_column: str):
    relevant_genes_and_tfs = [current_gene]
    tf_entries = row[name_of_jaspar_column]
    if pd.isna(tf_entries) or tf_entries == ".":
        return relevant_genes_and_tfs
    for entry in tf_entries.split(";"):
        tf_name = entry.split(",")[0]
        relevant_genes_and_tfs.append(tf_name)
    relevant_genes_and_tfs = list(set(relevant_genes_and_tfs))
    return relevant_genes_and_tfs


def _create_dictionary_for_pids_in_cohort(path_to_metadata):
    """
    cohort_dict[cohort1] = [pid1, pid2, ...]
    """
    metadata = pd.read_csv(path_to_metadata)
    primary_tumor_df = metadata[(metadata["cancer_type"] == "tumor")]
    primary_tumor_df.reset_index(inplace=True, drop=True)

    cohort_dict = primary_tumor_df.groupby(
        'cohort')['pid'].apply(list).to_dict()
    return cohort_dict


def _get_gene_and_tf_row_expressions(
    row: pd.Series,
    rnaseq_dataframe: pd.DataFrame,
    path_to_metadata: str,
    relevant_genes_and_tfs: list,
    relevant_cohorts: list
):
    # Get those genes/TFs that are within the relevant genes and tfs list and also within the RNA-seq dataframe.
    genes_and_tfs_in_rnaseq_dataframe = list(rnaseq_dataframe.columns)
    relevant_genes_and_tfs_in_rnaseq_dataframe = [
        gene_or_tf for gene_or_tf in relevant_genes_and_tfs if gene_or_tf in genes_and_tfs_in_rnaseq_dataframe]
    rnaseq_dataframe = rnaseq_dataframe[[
        gene_or_tf for gene_or_tf in relevant_genes_and_tfs if gene_or_tf in list(rnaseq_dataframe.columns)]]

    # Create a dictionary for each relevant cohort.
    # Ex: cohort_dict[cohort1] = [pid1, pid2, ...]
    cohort_dict = _create_dictionary_for_pids_in_cohort(path_to_metadata)

    # Create the row expressions dictionary.
    # Ex: row_expressions[tf_name][cohort][raw|zscore|log]
    row_expressions = {}
    for cohort in relevant_cohorts:
        relevant_pids = cohort_dict[cohort]
        raw_expression = rnaseq_dataframe.loc[[
            pid for pid in relevant_pids if pid in list(rnaseq_dataframe.index)]]
        if raw_expression.shape[0] == 0:
            continue
        zscore_expression = zscore(raw_expression, axis=0)
        log_expression = np.log(raw_expression)
        for gene_or_tf in relevant_genes_and_tfs_in_rnaseq_dataframe:
            if gene_or_tf not in row_expressions:
                row_expressions[gene_or_tf] = {}
            if cohort not in row_expressions[gene_or_tf]:
                row_expressions[gene_or_tf][cohort] = {
                    "raw": list(raw_expression[gene_or_tf]),
                    "zscore": list(zscore_expression[gene_or_tf]),
                    "log": list(log_expression[gene_or_tf])
                }

    return row_expressions


def _add_expression_traces_for_pSNV_hunter(
    dataframe: pd.DataFrame,
    pid: str,
    cohort_of_patient: str,
    recurrence_column_name: str,
    jaspar_column_name: str,
    path_to_metadata: str,
    dataset: str,
    path_to_rnaseq_dataframe: str,
):
    if "expression_traces" in list(dataframe.columns):
        print("Expression traces already exists in the dataframe!")
        return dataframe

    # Define the recurrence column.
    recurrence_column = [col for col in list(
        dataframe.columns) if col.endswith(recurrence_column_name)][0]

    # Read in the entire RNA-Seq dataframe and turn the pid into the index
    print("\tLoading in RNA-Seq dataframe...")
    rnaseq_dataframe = pd.read_csv(path_to_rnaseq_dataframe, delimiter="\t")
    rnaseq_dataframe.rename(columns={"Unnamed: 0": "pid"}, inplace=True)
    rnaseq_dataframe.set_index("pid", inplace=True)
    print("\t...done")

    # Iterate through each row in the dataframe.
    for idx, row in tqdm(dataframe.iterrows(), total=dataframe.shape[0], desc="Adding expression traces"):
        # Get the relevant cohorts (i.e. those cohorts that this SNV is recurrent with).
        relevant_cohorts = _get_relevant_cohorts_for_single_row(
            row, cohort_of_patient, recurrence_column)

        # Get the relevant genes/TFs (i.e. those genes/TFs that are affected by this SNV).
        relevant_genes_and_tfs = _get_relevant_genes_and_tfs_for_single_row(
            row, row["GENE"], jaspar_column_name)

        # Get the row expressions.
        # Format: row_expressions[tf_name][cohort][raw|zscore|log]
        row_expressions = _get_gene_and_tf_row_expressions(
            row=row,
            rnaseq_dataframe=rnaseq_dataframe,
            path_to_metadata=path_to_metadata,
            relevant_genes_and_tfs=relevant_genes_and_tfs,
            relevant_cohorts=relevant_cohorts
        )

        # Parse this dictionary into a string for writeability.
        dataframe.loc[idx, "expression_traces"] = str(
            row_expressions).replace("inf", "9999")

    return dataframe


def _sum_the_number_of_recurrent_mutations(
    dataframe: pd.DataFrame,
):
    recurrence_columns = [col for col in list(
        dataframe.columns) if col.endswith("num_recurrent_mutations")]
    dataframe["total_num_recurrent_mutations"] = np.sum(
        dataframe[recurrence_columns], axis=1)
    return dataframe


def _add_remind_cancer_score(
    dataframe: pd.DataFrame,
    remind_cancer_scoring_weights: dict,
    ge_column: str
):
    for idx, row in tqdm(dataframe.iterrows(), total=dataframe.shape[0], desc="Adding REMIND-Cancer Score"):
        # (1) Genomic.
        # Transcription Factors
        created_tf_weight = np.min([
            float(row["num_created_tfs_passing_tf_expression_threshold"]) *
            remind_cancer_scoring_weights["genomic"]["tfbs"]["creation_weight_per_tfbs"],
            remind_cancer_scoring_weights["genomic"]["tfbs"]["creation_weight_maximum"]
        ])

        destroyed_tf_weight = np.min([
            float(row["num_destroyed_tfs_passing_tf_expression_threshold"]) *
            remind_cancer_scoring_weights["genomic"]["tfbs"]["destruction_weight_per_tfbs"],
            remind_cancer_scoring_weights["genomic"]["tfbs"]["destruction_weight_maximum"]
        ])

        # Recurrence
        recurrence_weight = np.min([
            float(row["total_num_recurrent_mutations"]) *
            remind_cancer_scoring_weights["genomic"]["recurrence"]["weight_per_recurrent_mutation"],
            remind_cancer_scoring_weights["genomic"]["recurrence"]["weight_maximum"]
        ])

        # Purity
        try:
            purity_weight = remind_cancer_scoring_weights["genomic"]["purity"]["purity_weight"] if float(
                row["purity"]) >= remind_cancer_scoring_weights["genomic"]["purity"]["purity_threshold_for_weight"] else 0
        except:
            print(
                f"Purity of {row['purity']} does not work. Setting score to 0.")
            purity_weight = 0

        # Allele Frequency
        try:
            af_weight = remind_cancer_scoring_weights["genomic"]["allele_frequency"]["af_weight"] if float(
                row["allele_frequency"]) >= remind_cancer_scoring_weights["genomic"]["allele_frequency"]["af_threshold_for_weight"] else 0
        except:
            print(
                f"Allele frequency of {row['allele_frequency']} does not work. Setting score to 0.")
            af_weight = 0

        # (2) Transcriptomic
        # Gene Expression
        gene_expression_weight = np.min([
            float(row[ge_column]) *
            remind_cancer_scoring_weights["transcriptomic"]["gene_expression"]["weight_per_unit_of_expression"],
            remind_cancer_scoring_weights["transcriptomic"]["gene_expression"]["weight_maximum"]
        ])

        # (3) Annotations
        # CGC
        cgc_weight = bool(row["within_cgc_list"]) * \
            remind_cancer_scoring_weights["annotations"]["cgc"]["weight"]

        # Open Chromatin
        open_chromatin_weight = bool(
            row["open_chromatin"]) * remind_cancer_scoring_weights["annotations"]["open_chromatin"]["weight"]

        dataframe.loc[idx, "score"] = np.sum([
            created_tf_weight,
            destroyed_tf_weight,
            recurrence_weight,
            purity_weight,
            af_weight,
            gene_expression_weight,
            cgc_weight,
            open_chromatin_weight
        ])

    dataframe = dataframe.sort_values("score", ascending=False)

    return dataframe


def main(
    path_to_vcf_file: str,
    config: dict,
):
    # Read in the dataframe.
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")

    # Add the PID to the dataframe.
    pid = _get_pid_from_structured_vcf_path(path_to_vcf_file, True)
    data["pid"] = pid

    # Add the cohort to the dataframe.
    patient_to_cohort_dictionary = _create_cohort_dictionary(
        config["pipeline"]["path_to_metadata"])
    cohort_of_patient = patient_to_cohort_dictionary[pid]
    data["cohort"] = cohort_of_patient

    # Add the TF logo plot.
    print("Adding TF logo plots...")
    path_to_downloaded_jaspar_json = config["post_processing"][
        "transcription_factor_information"]["path_to_downloaded_jaspar_file"]
    data = _add_tf_logo_sequence_to_single_patient_path(
        dataframe=data,
        path_to_downloaded_jaspar_json=path_to_downloaded_jaspar_json,
    )

    # Add the number of original, promoter and final mutations.
    print("Adding the number of original, promoter, and final mutations.")
    with open(config["pipeline"]["path_to_results"], "r") as f:
        results = json.load(f)

    ending_after_preprocessing = config["preprocessing_details"]["suffix_to_append_to_vcf"]
    ending_after_promoter = config["promoter"]["suffix_to_append_to_vcf"]

    data = _add_num_original_promoter_and_final_mutations_to_single_patient_path(
        dataframe=data,
        pid=pid,
        results=results,
        ending_after_preprocessing=ending_after_preprocessing,
        ending_after_promoter=ending_after_promoter
    )

    # Add NCBI information.
    print("Adding NCBI information to genes and transcrciption factors.")
    path_to_gene_name_and_description_file = config["post_processing"][
        "ncbi_information"]["path_to_gene_name_and_description_file"]
    name_of_tfbs_prediction_column = config["preprocessing_details"][
        "transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"]

    data = _add_ncbi_information_of_genes_and_tfs(
        dataframe=data,
        path_to_gene_name_and_description_file=path_to_gene_name_and_description_file,
        name_of_tfbs_prediction_column=name_of_tfbs_prediction_column
    )

    # Add expression traces for pSNV hunter.
    print("Adding expression traces for pSNV hunter.")
    recurrence_column_name = config["recurrence"]["name_of_column_to_add"]
    jaspar_column_name = config["preprocessing_details"]["transcription_factor_prediction"][
        "name_of_tfbs_prediction_column_to_add"] + "(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"

    path_to_metadata = config["pipeline"]["path_to_metadata"]
    dataset = config["pipeline"]["dataset"]
    path_to_rnaseq_dataframe = config["rnaseq_data"][dataset]["path_to_rnaseq_dataframe"]

    patient_cohort_dictionary = _create_cohort_dictionary(path_to_metadata)

    # data = _add_expression_traces_for_pSNV_hunter(
    #     dataframe=data,
    #     pid=pid,
    #     cohort_of_patient=cohort_of_patient,
    #     recurrence_column_name=recurrence_column_name,
    #     jaspar_column_name=jaspar_column_name,
    #     path_to_metadata=path_to_metadata,
    #     dataset=dataset,
    #     path_to_rnaseq_dataframe=path_to_rnaseq_dataframe,
    # )

    # Add the total number of recurrent mutations.
    print("Adding the _total_ number of recurrent mutations.")
    data = _sum_the_number_of_recurrent_mutations(data)

    # Add the REMIND-Cancer score.
    print("Adding the REMIND-Cancer score.")
    ge_column = config["ge_filter"]["column_name_to_filter"]
    remind_cancer_scoring_weights = config["REMIND-Cancer_scoring_weights"]

    data = _add_remind_cancer_score(
        dataframe=data,
        remind_cancer_scoring_weights=remind_cancer_scoring_weights,
        ge_column=ge_column
    )

    # Write the file.
    print("Writing the new file.")
    suffix_to_append_to_vcf = config["post_processing"]["suffix_to_append_to_vcf"]
    data.to_csv(
        path_to_vcf_file.replace(".vcf", suffix_to_append_to_vcf),
        sep="\t", index=False
    )

    print(
        f"New postprocessing file written to {path_to_vcf_file.replace('.vcf', suffix_to_append_to_vcf)}!")


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
        name_of_pipeline_step = "postprocessing",
        start_time = start_time,
        end_time = end_time,
        name_of_timing_tracker_file = "timing_tracker.json",
    )
    
    print("~~Finished with postprocessing~~")
