import pandas as pd
import os
from optparse import OptionParser
import json

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)
    
from src.pipeline.general_helper_functions import execution_time



@execution_time
def _ge_filter(
    path_to_vcf_file: str,
    threshold: float,
    column_name_to_filter: str,
    suffix_to_append_to_vcf: str,
):
    # Load in the files.
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")

    # Remove those that do not have a corresponding Z-score.
    data = data[data[column_name_to_filter] != "not_available"]

    # If there are any rows with the Z-score being NaN, set to 0.
    data[column_name_to_filter] = data[column_name_to_filter].fillna(0)

    # If there are no remaining rows, return nothing.
    if data.shape[0] == 0:
        return

    # Reset the index.
    data.reset_index(drop=True, inplace=True)

    # Remove those that are not above the threshold.
    data = data[data[column_name_to_filter] > threshold]
    if data.shape[0] == 0:
        return

    # Reset the index.
    data.reset_index(drop=True, inplace=True)

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
    threshold = config["ge_filter"]["threshold"]
    column_name_to_filter = config["ge_filter"]["column_name_to_filter"]
    suffix_to_append_to_vcf = config["ge_filter"]["suffix_to_append_to_vcf"]

    _ge_filter(
        path_to_vcf_file=path_to_vcf_file,
        threshold=threshold,
        column_name_to_filter=column_name_to_filter,
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

    main(
        path_to_vcf_file=path_to_vcf_file,
        config=config
    )

    print("~~ Finished with GE filter ~~")
