import pandas as pd
from optparse import OptionParser
import json

import sys
import os

# For proper import structure:
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import execution_time


@execution_time
def _add_cgc_information(
    path_to_vcf_file: str,
    path_to_cgc_file: str,
    suffix_to_append_to_vcf: str,
):
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")
    cgc = pd.read_csv(path_to_cgc_file)

    cgc_genes = list(cgc["Gene Symbol"].unique())

    data["within_cgc_list"] = data["GENE"].isin(cgc_genes)

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
    path_to_cgc_file = config["cgc"]["path_to_cgc_file"]
    suffix_to_append_to_vcf = config["cgc"]["suffix_to_append_to_vcf"]

    _add_cgc_information(
        path_to_vcf_file=path_to_vcf_file,
        path_to_cgc_file=path_to_cgc_file,
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

    print("~~ Finished with adding CGC information ~~")
