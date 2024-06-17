import pandas as pd
from tqdm import tqdm
from optparse import OptionParser
import json
import os

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import execution_time

@execution_time
def _promoter_filter(
    path_to_vcf_file: str,
    path_to_tss_file: str,
    upstream_bp: int,
    downstream_bp: int,
    suffix_to_append_to_vcf: str,
):
    # Load in the files.
    tss_file = pd.read_csv(path_to_tss_file, delimiter="\t")
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")

    # Create the promoter start and end.
    for idx, row in tqdm(tss_file.iterrows(), total=tss_file.shape[0]):
        strand = row["strand"]
        if strand == "+":
            tss_file.loc[idx, "promoter_start"] = row["gene_start"] - upstream_bp
            tss_file.loc[idx, "promoter_end"] = row["gene_start"] + downstream_bp
        elif strand == "-":
            tss_file.loc[idx, "promoter_start"] = row["gene_end"] - downstream_bp
            tss_file.loc[idx, "promoter_end"] = row["gene_end"] + upstream_bp 

    # Check whether a row (1) is in the chromosome, (2) same gene name, (3) in between the promoter start and end.
    for idx, row in tqdm(data.iterrows(), total=data.shape[0]):
        vcf_chromosome = f"chr{row['#CHROM']}"
        vcf_position = int(row["POS"])
        vcf_gene = str(row["GENE"])

        within_promoter = tss_file[
            (tss_file["#chr"] == vcf_chromosome) & 
            (tss_file["gene_name"] == vcf_gene) &
            (tss_file["promoter_start"] <= vcf_position) &
            (tss_file["promoter_end"] >= vcf_position)
        ].shape[0] > 0

        data.loc[idx, "within_promoter"] = within_promoter

    # Turn data into those that are within the promoter.
    data = data[data["within_promoter"] == True]
    
    if data.shape[0] == 0:
        return
    
    # Drop the "within_promoter" column
    data.drop("within_promoter", axis=1, inplace=True)

    data.reset_index(
        inplace=True,
        drop=True,
    )
    
    data.to_csv(
        path_to_vcf_file.replace(".vcf", suffix_to_append_to_vcf),
        sep="\t",
        index=False
    )

@execution_time
def main(
    path_to_vcf_file: str,
    config: dict
):
    path_to_tss_file = config["additional_files"]["tss_reference_file"]
    upstream_bp = config["promoter"]["upstream_of_tss"]
    downstream_bp = config["promoter"]["downstream_of_tss"]
    suffix_to_append_to_vcf = config["promoter"]["suffix_to_append_to_vcf"]
    
    _promoter_filter(
        path_to_vcf_file = path_to_vcf_file,
        path_to_tss_file = path_to_tss_file,
        upstream_bp = upstream_bp,
        downstream_bp = downstream_bp,
        suffix_to_append_to_vcf = suffix_to_append_to_vcf,
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
        path_to_vcf_file = path_to_vcf_file,
        config = config
    )
    
    print("~~ Finished with promoter filter ~~")
    