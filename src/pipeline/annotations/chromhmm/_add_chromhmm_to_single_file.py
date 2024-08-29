import pandas as pd
from tqdm import tqdm
from optparse import OptionParser
import json
from datetime import datetime

# For proper import structure:
import os
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)
    
from src.pipeline.general_helper_functions import execution_time, _update_timing_tracker

# Function to check if position falls within any of the ranges for a chromosome
def check_if_snv_within_range(row: pd.Series, chromhmm_df: pd.DataFrame):
    chromosome = str(row['#CHROM'])
    position = int(row['POS'])

    relevant_rows = chromhmm_df[chromhmm_df['chromosome'].astype(
        str) == chromosome]

    for _, chrom_row in relevant_rows.iterrows():
        if chrom_row['start_pos'] <= position <= chrom_row['end_pos']:
            return True
    return False


@execution_time
def _add_chromatin_information(
    path_to_vcf_file: str,
    path_to_chromhmm_file: str,
    suffix_to_append_to_vcf: str,
):

    data = pd.read_csv(path_to_vcf_file, delimiter="\t")
    chromhmm = pd.read_csv(path_to_chromhmm_file, delimiter="\t", names=[
                           'chromosome', 'start_pos', 'end_pos', 'classification', 'classification2'])

    # Ensure that both chromosome columns are in the same format.
    chromhmm["chromosome"] = chromhmm["chromosome"].str.replace(
        "chr", "").astype(str)
    data["#CHROM"] = data["#CHROM"].astype(str)

    # Apply the check_if_snv_within_range function to the dataset.
    tqdm.pandas(desc="Checking ranges")
    data['open_chromatin'] = data.progress_apply(
        check_if_snv_within_range, axis=1, args=(chromhmm,))

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
    suffix_to_append_to_vcf = config["open_chromatin"]["suffix_to_append_to_vcf"]
    path_to_chromhmm_file = config["open_chromatin"]["path_to_chromhmm_file"]

    _add_chromatin_information(
        path_to_vcf_file=path_to_vcf_file,
        path_to_chromhmm_file=path_to_chromhmm_file,
        suffix_to_append_to_vcf=suffix_to_append_to_vcf,
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
        name_of_pipeline_step = "chromhmm",
        start_time = start_time,
        end_time = end_time,
        name_of_timing_tracker_file = "timing_tracker.json",
    )

    print("~~ Finished with adding open chromatin information ~~")
