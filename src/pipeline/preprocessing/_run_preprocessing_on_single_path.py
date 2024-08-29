import pandas as pd
import os
from tqdm import tqdm
import shutil
import re
import numpy as np
import glob

from datetime import datetime
import json
import subprocess

from optparse import OptionParser
from scipy.stats import zscore
from optparse import OptionParser
from pyfaidx import Fasta

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import _get_index, _get_pid_from_structured_vcf_path, execution_time, _update_timing_tracker

@execution_time
def _pcawg_rename_hugo_symbol_to_gene(
    path_to_single_vcf: str
):
    data = pd.read_csv(path_to_single_vcf, delimiter="\t")

    if "Hugo_Symbol" in list(data.columns):
        data.rename({"Hugo_Symbol": "GENE"}, axis=1, inplace=True)
        data.to_csv(path_to_single_vcf, sep="\t", index=False)

    return path_to_single_vcf


@execution_time
def _fix_problematic_genes(path_to_single_vcf: str):
    """Fix the problem in regards to genes and their names.

    To exemplify this problem, in the original .vcf files, there are mutations within genes that have the following names:
    (1) "EMILIN1(ENST00000380320.4:c.-28T>A)"
    (2) "TNPO1(ENST00000337273.5:c.-124C>G,ENST00000454282.1:c.-124C>G)"
    (3) "NDUFB11(ENST00000377811.3:c.-753G>T),RBM10(ENST00000377604.3:c.-2050C>A,ENST00000329236.7:c.-2050C>A,ENST00000345781.6:c.-2050C>A)"
    (4) "COX6A1P2(dist=80875),PIM1(dist=43927)"

    The information within the parenthesis are stripped and if there are multiple genes, a line for each gene is created with a corresponding gene.

    For example, here would be the gene names of the above examples:
    (1) "EMILIN1"
    (2) "TNPO1"
    (3) "NDUFB11,RBM10" -> (results in 2 separate lines)
    (4) "COX6A1P2,PIM1" -> (results in 2 separate lines)

    """ ""
    assert os.path.exists(
        path_to_single_vcf), f"_fix_problematic_genes -> Cannot read in {path_to_single_vcf}"

    data = pd.read_csv(path_to_single_vcf, delimiter="\t")

    # Remove everything within parenthesis.
    # Ex: CCDC132(dist=50563),CALCR(dist=14898) -> CCDC132,CALCR
    # Ex: IQCE(ENST00000404984.1:c.*2600A>G,ENST00000402050.2:c.*2600A>G) -> IQCE
    data["GENE"] = data["GENE"].apply(lambda x: re.sub(r"\([^)]*\)", "", x))

    # In the case that there are multiple genes, replicate the line for each of these genes.
    # Ex: "NDUFB11,RBM10" goes to two lines with the same exact info except line 1 has NDUFB11 as its gene and line 2 has RBM10
    expanded_gene_column = data["GENE"].str.split(
        ',').explode().str.split(';').explode()

    # Concat them along the index.
    data = pd.concat(
        [data.drop("GENE", inplace=False, axis=1), expanded_gene_column], axis=1)

    # Reset the index.
    data.reset_index(drop=True, inplace=True)

    # Re-write to the original location.
    data.to_csv(
        path_to_single_vcf,
        sep="\t", index=False
    )

    return path_to_single_vcf


@execution_time
def _remove_lines_with_no_gene(path_to_single_vcf: str):
    # Read in the dataframe.
    data = pd.read_csv(path_to_single_vcf, delimiter="\t")

    # Remove those genes that are labeled as "NONE".
    # "NONE" could happen in MASTER, "Unknown" can happen in PCAWG.
    data = data[(data["GENE"] != "NONE") & (data["GENE"] != "Unknown")]

    # Remove those that are NaN in the GENE column.
    data["GENE"].dropna(inplace=True)

    # Reset the index.
    data.reset_index(drop=True, inplace=True)

    data.to_csv(
        path_to_single_vcf,
        sep="\t", index=False
    )

    return path_to_single_vcf


@execution_time
def _remove_indels(path_to_single_vcf: str):

    # Read in the dataframe.
    data = pd.read_csv(path_to_single_vcf, delimiter="\t")

    # Remove those such that the REF or ALT do not have a length of 1.
    data = data[(data['REF'].str.len() == 1) & (data['ALT'].str.len() == 1)]

    # Reset index.
    data.reset_index(drop=True, inplace=True)

    data.to_csv(path_to_single_vcf, sep="\t", index=False)

    return path_to_single_vcf


@execution_time
def _remove_snp_and_germline(path_to_single_vcf: str):

    # Read in the dataframe.
    data = pd.read_csv(path_to_single_vcf, delimiter="\t")

    # Within DKFZ's SNV Calling Workflow (https://github.com/DKFZ-ODCF/SNVCallingWorkflow/), there's a column called
    # "RECLASSIFICATION". If labeled "somatic", they're somatic whereas the others are likely not.

    if "RECLASSIFICATION" in list(data.columns):
        data = data[data["RECLASSIFICATION"] == "somatic"]
        data.reset_index(drop=True, inplace=True)
        data.to_csv(path_to_single_vcf, index=False, sep="\t")

    return path_to_single_vcf


@execution_time
def _add_sequence_context_column_to_pcawg_vcf_files(
    path: str
):
    """Because the original PCAWG data files do not have the SEQUENCE_CONTEXT column, create this manually.

    To do so, we take the 'ref_context' column for each line and get the nucleotide that is directly in the middle.
    With this nucleotide, we make sure that this is equal to the reference nucleotide.

    With these, we then re-create the SEQUENCE_CONTEXT column by replacing the ref_context at the middle position
    with a "," to emulate the same structure as the MASTER dataset.
    """
    with open(path, "r") as input:
        original_lines = input.readlines()

    if "SEQUENCE_CONTEXT" in original_lines[0]:
        return path

    with open(path, "w") as output:
        # Write the new header to the output file.
        input_header = original_lines[0]
        output_header = input_header.replace("\n", "\t") + "SEQUENCE_CONTEXT\n"
        output.write(output_header)

        # Get the necessary indices.
        ref_idx = _get_index(input_header, "REF")
        ref_context_idx = _get_index(input_header, "ref_context")

        for line in original_lines[1:]:
            line_split = line.split("\t")
            ref = line_split[ref_idx].upper()
            ref_context = line_split[ref_context_idx].upper()
            position_of_ref = int(len(ref_context) / 2)
            if ref == ref_context[position_of_ref]:
                sequence_context = ""
                for idx, nucleotide in enumerate(ref_context):
                    if idx == position_of_ref:
                        sequence_context += ","
                    else:
                        sequence_context += nucleotide
                output.write(
                    line.replace("\n", "\t") + sequence_context + "\n"
                )
    return path


def _get_fimo_output(
    path_to_vcf_file: str,
    path_to_fimo: str,
    path_to_jaspar_database: str
):
    # Create temporary fasta files, which will be deleted at the end of the script.
    # These will be in the same folder as the original .vcf file.
    temporary_ref_fasta_file = path_to_vcf_file.replace(
        ".vcf", "_temp_ref_fasta.fna"
    )

    temporary_alt_fasta_file = path_to_vcf_file.replace(
        ".vcf", "_temp_alt_fasta.fna"
    )

    # Define the output of FIMO.
    ref_output_file = path_to_vcf_file.replace(
        ".vcf", "_ref_output.tsv"
    )
    alt_output_file = path_to_vcf_file.replace(
        ".vcf", "_alt_output.tsv"
    )

    # Using the inputted .vcf, write two fasta files (ref and alt) that will be used for FIMO.
    with open(path_to_vcf_file, "r") as f, open(temporary_ref_fasta_file, "w") as ref_fasta_file, open(temporary_alt_fasta_file, "w") as alt_fasta_file:
        # Read in the header.
        header = f.readline()

        # Get the indices of the important columns.
        chrom_idx = _get_index(header, "#CHROM", True)
        pos_idx = _get_index(header, "POS", True)
        ref_idx = _get_index(header, "REF", True)
        alt_idx = _get_index(header, "ALT", True)
        sequence_context_idx = _get_index(header, "SEQUENCE_CONTEXT", True)

        # Iterate through each line in the input dataframe and put this into fasta format.
        # This is done for both the ref and alt.
        for idx, line in enumerate(f):
            line_split = line.split("\t")
            chrom = line_split[chrom_idx].strip()
            pos = line_split[pos_idx].strip()
            ref = line_split[ref_idx].strip()
            alt = line_split[alt_idx].strip()
            sequence_context = line_split[sequence_context_idx].strip()

            # Write this line to both the reference and the alt fasta files.
            ref_fasta_file.write(f">{idx}_ref_{chrom}_{pos}_{ref}_{alt}\n")
            ref_fasta_file.write(f"{sequence_context.replace(',', ref)}\n")

            alt_fasta_file.write(f">{idx}_alt_{chrom}_{pos}_{ref}_{alt}\n")
            alt_fasta_file.write(f"{sequence_context.replace(',', alt)}\n")

    # Call FIMO on both of the files.
    print(f"Working on ref...")
    command = f"{path_to_fimo} --verbosity 1 --text {path_to_jaspar_database} {temporary_ref_fasta_file} > {ref_output_file}"
    process = subprocess.Popen(command, shell=True)
    process.wait()

    print(f"Working on alt...")
    command = f"{path_to_fimo} --verbosity 1 --text {path_to_jaspar_database} {temporary_alt_fasta_file} > {alt_output_file}"
    process = subprocess.Popen(command, shell=True)
    process.wait()

    # Try to remove the temporary fasta files.
    try:
        os.remove(temporary_ref_fasta_file)
        os.remove(temporary_alt_fasta_file)
    except:
        pass

    return ref_output_file, alt_output_file


def _process_output_from_fimo(
    ref_output_file: str,
    alt_output_file: str,
):
    # Process the output files computed in _get_fimo_output().
    with open(ref_output_file, "r") as ref_output, open(alt_output_file, "r") as alt_output:
        # Get the headers, which should be the same.
        ref_header = ref_output.readline()
        alt_header = alt_output.readline()

        # Get the necessary column indices from the output of FIMO.
        tf_name_idx = _get_index(ref_header, "motif_alt_id", True)
        motif_id_idx = _get_index(ref_header, "motif_id", True)
        sequence_name_idx = _get_index(ref_header, "sequence_name", True)
        score_idx = _get_index(ref_header, "score", True)
        matched_sequence_idx = _get_index(ref_header, "matched_sequence", True)

        # Create a dictionary to keep track of TF predictions.
        # Ex: tfbs_dictionary[0]["ELK4"]["ref"]["score"] = 10
        # Ex: tfbs_dictionary[0]["ELK4"]["alt"]["score"] = 5
        # Ex: tfbs_dictionary[0]["ELK4"]["ratio"] = 0.5
        tfbs_dictionary = {}

        # Add the reference file lines to the dictionary.
        for line in tqdm(ref_output, desc="Adding ref to dictionary..."):
            line_split = line.split("\t")

            try:
                tf_name = line_split[tf_name_idx].strip(
                )                   # Ex: Alx1
                # Ex: MA0854.1
                motif_id = line_split[motif_id_idx].strip()
                # Ex: 1231_alt_3_183393144_G_T
                sequence_name = line_split[sequence_name_idx].strip()
                # Ex: 4.33e-05
                score = float(str(line_split[score_idx].strip()))
                # Ex: CCCTTTAATTAACTGGG
                matched_sequence = line_split[matched_sequence_idx].strip()
            except:
                print(f"REMIND-Cancer: Error with {line_split}")
                continue

            corresponding_index_in_vcf_file = int(sequence_name.split("_")[0])

            if corresponding_index_in_vcf_file not in tfbs_dictionary:
                tfbs_dictionary[corresponding_index_in_vcf_file] = {}

            if tf_name not in tfbs_dictionary[corresponding_index_in_vcf_file]:
                tfbs_dictionary[corresponding_index_in_vcf_file][tf_name] = {
                    "alt": {"score": 1, "sequence": "."},
                    "ref": {"score": 1, "sequence": "."}
                }

            tfbs_dictionary[corresponding_index_in_vcf_file][tf_name]["ref"]["score"] = score
            tfbs_dictionary[corresponding_index_in_vcf_file][tf_name]["ref"]["sequence"] = matched_sequence

        # Add the alt file lines to the dictionary.
        for line in tqdm(alt_output, desc="Adding alt to dictionary..."):
            line_split = line.split("\t")

            try:
                tf_name = line_split[tf_name_idx].strip(
                )                   # Ex: Alx1
                # Ex: MA0854.1
                motif_id = line_split[motif_id_idx].strip()
                # Ex: 1231_alt_3_183393144_G_T
                sequence_name = line_split[sequence_name_idx].strip()
                # Ex: 4.33e-05
                score = float(str(line_split[score_idx].strip()))
                # Ex: CCCTTTAATTAACTGGG
                matched_sequence = line_split[matched_sequence_idx].strip()
            except:
                print(f"REMIND-Cancer: Error with {line_split}")
                continue

            corresponding_index_in_vcf_file = int(sequence_name.split("_")[0])

            if corresponding_index_in_vcf_file not in tfbs_dictionary:
                tfbs_dictionary[corresponding_index_in_vcf_file] = {}

            if tf_name not in tfbs_dictionary[corresponding_index_in_vcf_file]:
                tfbs_dictionary[corresponding_index_in_vcf_file][tf_name] = {
                    "alt": {"score": 1, "sequence": "."},
                    "ref": {"score": 1, "sequence": "."}
                }

            tfbs_dictionary[corresponding_index_in_vcf_file][tf_name]["alt"]["score"] = score
            tfbs_dictionary[corresponding_index_in_vcf_file][tf_name]["alt"]["sequence"] = matched_sequence

        # Calculate the scoring ratio between the reference and the alts. (alt / ref)
        # Also, get the entry of this TF in the proper format. (tf_name,binding_affinity,seq1,seq2)
        for idx in tfbs_dictionary:
            for tf_name in tfbs_dictionary[idx]:
                try:
                    tfbs_dictionary[idx][tf_name]["ratio"] = np.abs(
                        tfbs_dictionary[idx][tf_name]["alt"]["score"] /
                        tfbs_dictionary[idx][tf_name]["ref"]["score"]
                    )
                except:
                    tfbs_dictionary[idx][tf_name]["ratio"] = 0
                tfbs_dictionary[idx][tf_name]["string"] = f"{tf_name},{tfbs_dictionary[idx][tf_name]['ratio']},{tfbs_dictionary[idx][tf_name]['ref']['sequence']},{tfbs_dictionary[idx][tf_name]['alt']['sequence']}"
    return tfbs_dictionary


def _add_processed_fimo_output_to_vcf(
    path_to_vcf_file: str,
    tfbs_dictionary: dict,
    name_of_jaspar_column: str,
):
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")

    for idx in tfbs_dictionary:
        jaspar_entry = ""
        for tf_name in tfbs_dictionary[idx]:
            jaspar_entry += tfbs_dictionary[idx][tf_name]["string"] + ";"
        jaspar_entry = jaspar_entry[:-1]
        data.loc[idx, name_of_jaspar_column] = jaspar_entry

    data.fillna({name_of_jaspar_column: "."}, inplace=True)
    return data


@execution_time
def run_fimo(
    path_to_vcf_file: str,
    path_to_fimo: str = "/home/n795d/meme/bin/fimo",
    path_to_jaspar_database: str = "/home/n795d/workspace/REMIND-Cancer/data/general_input_files/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt",
    name_of_jaspar_column: str = "JASPAR2020_CORE_vertebrates_non-redundant",
    save_tfbs_dict: bool = True,
):
    assert os.path.exists(
        path_to_vcf_file), "Path to .vcf file does not exist."
    assert os.path.exists(path_to_fimo), "Path to fimo file does not exist."
    assert os.path.exists(
        path_to_jaspar_database), "Path to jaspar database file does not exist."

    ref_output_file, alt_output_file = _get_fimo_output(
        path_to_vcf_file=path_to_vcf_file,
        path_to_fimo=path_to_fimo,
        path_to_jaspar_database=path_to_jaspar_database
    )

    tfbs_dictionary = _process_output_from_fimo(
        ref_output_file=ref_output_file,
        alt_output_file=alt_output_file
    )

    data = _add_processed_fimo_output_to_vcf(
        path_to_vcf_file=path_to_vcf_file,
        tfbs_dictionary=tfbs_dictionary,
        name_of_jaspar_column=name_of_jaspar_column
    )

    if save_tfbs_dict:
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime("%d-%m-%Y_%H:%M")

        output_path_of_saved_tfbs_dict = path_to_vcf_file.replace(
            ".vcf", f"_fimo_dict_{formatted_datetime}.json")

        print(f"Saving TFBS dict to {output_path_of_saved_tfbs_dict}")
        with open(output_path_of_saved_tfbs_dict, "w") as f:
            json.dump(tfbs_dictionary, f)
        print(f"... done.")
    data.to_csv(
        path_to_vcf_file, sep="\t", index=False
    )

    return path_to_vcf_file


@execution_time
def _fix_transcription_factor_names(path_to_single_vcf: str, name_of_jaspar_column: str = "JASPAR2020_CORE_vertebrates_non-redundant"):
    with open(path_to_single_vcf, "r") as f:
        header = f.readline()
        jaspar_idx = _get_index(header, name_of_jaspar_column, True)
        new_lines = f"{header}"
        for line in f:
            line_split = line.split("\t")
            jaspar = line_split[jaspar_idx]
            fixed_entry = ""
            for entry in jaspar.split(";"):
                try:
                    gene, binding_affinity, extra_info1, extra_info2 = entry.split(
                        ",")
                    binding_affinity = float(str(binding_affinity).strip())
                    extra_info1 = extra_info1.strip()
                    extra_info2 = extra_info2.strip()

                except:
                    continue

                # Remove parenthesis within TF name (i.e. JUND(var.2) -> JUND)
                gene = re.sub(r"\([^)]*\)", "", gene)

                # Check if there are any lowercase letters.
                lower_case_letter_exists = any(
                    letter.islower() for letter in gene)
                if lower_case_letter_exists:
                    continue

                # Split genes that are two. (i.e. FOSL1::JUND -> FOSL1 and JUND)
                if "::" in gene:
                    gene1 = gene.split("::")[0]
                    gene2 = gene.split("::")[1]

                    fixed_entry += f"{gene1},{binding_affinity},{extra_info1},{extra_info2};{gene2},{binding_affinity},{extra_info1},{extra_info2};"

                else:
                    fixed_entry += f"{gene},{binding_affinity},{extra_info1},{extra_info2};"
            if fixed_entry == "":
                fixed_entry = ".;"

            fixed_entry = fixed_entry[:-1]

            line_split[jaspar_idx] = fixed_entry

            new_line = "\t".join(line_split)

            if not new_line.endswith("\n"):
                new_line = new_line + "\n"

            new_lines += new_line

    with open(path_to_single_vcf, "w") as f:
        f.write(new_lines)

    return path_to_single_vcf


@execution_time
def _add_vaf_information_to_final_vcf_files(path_to_vcf, dataset):
    vcf_data = pd.read_csv(path_to_vcf, delimiter="\t")
    if "allele_frequency" in vcf_data.columns:
        return
    if dataset == "master":
        for idx, row in vcf_data.iterrows():
            dp4_info = row["INFO"].split("DP4=")[1].split(";")[0].split(",")
            dp4_info = [float(i) for i in dp4_info]
            ref_total = dp4_info[0] + dp4_info[1]
            alt_total = dp4_info[2] + dp4_info[3]
            allele_freq = np.round(
                alt_total / (ref_total + alt_total), decimals=5
            )
            vcf_data.loc[idx, "allele_frequency"] = allele_freq
    elif dataset == "pcawg":
        vcf_data["allele_frequency"] = vcf_data["i_VAF"]

    vcf_data.to_csv(path_to_vcf, index=None, sep="\t")
    return path_to_vcf


@execution_time
def _add_purity_information(
    path_to_vcf: str,
    config: dict
):
    dataset = config["pipeline"]["dataset"]

    if dataset == "master":
        path_to_metadata = config["pipeline"]["path_to_metadata"]

        # (1) Get the directory where the original SNV .vcf file comes from.
        # Within this directory is where the necessary plots are located.
        metadata = pd.read_csv(path_to_metadata)

        pid, cancer_type = _get_pid_from_structured_vcf_path(
            path_to_vcf).split("_")

        metadata_row_of_interest = metadata[(metadata["pid"] == pid) & (
            metadata["cancer_type"] == cancer_type)]

        assert metadata_row_of_interest.shape[0] == 1

        directory_of_original_vcf_file = os.path.dirname(
            metadata_row_of_interest["path_to_vcf"].iloc[0])

        # (2) Get the purity information if it exists.
        purity_file = None
        for file in glob.glob(f"{directory_of_original_vcf_file}/*"):
            # If the purity file exists, open the file and get the scores.
            if file.endswith("purityEST.txt"):
                purity_file = file
                median_after_esd, sd_after_esd, median_before_esd = _parse_master_file_for_purity(
                    file
                )

        if not purity_file:
            median_after_esd = "not_available"
            purity_file = "no_purity_file_exists"

        # (3) Open the original file and read the lines.
        with open(path_to_vcf, "r") as original:
            original_lines = original.readlines()
            # Check if the header already has the purity information. If so, skip this.
            if ("purity" in original_lines[0]):
                return path_to_vcf

        # (4) Re-write the original file with the added information.
        with open(path_to_vcf, "w") as new:
            new_header = original_lines[0].replace("\n", "\t")
            new_header += "purity\tpath_to_purity_file\n"
            new.write(new_header)

            for line in original_lines[1:]:
                new.write(
                    line.replace("\n", "\t")
                    + f"{median_after_esd}\t{purity_file}\n"
                )

    elif dataset == "pcawg":
        # Define the purity files.
        path_to_pcawg_purity_file = config["additional_files"]["pcawg"]["path_to_purity_file"]
        assert os.path.exists(
            path_to_pcawg_purity_file), "PCAWG Purity file does not exist."

        # Get the PCAWG pid of this file.
        pid = _get_pid_from_structured_vcf_path(path_to_vcf, True)

        # Load in the dataset and the purity file.
        data = pd.read_csv(path_to_vcf, delimiter="\t")
        purity_df = pd.read_csv(path_to_pcawg_purity_file, delimiter="\t")

        # Check if purity already exists as a column.
        if "purity" in list(data.columns):
            return path_to_vcf

        # Within the purity dataframe, find the corresponding pid if it exists.
        if pid in list(purity_df["samplename"]):
            purity = purity_df[purity_df["samplename"]
                               == pid].iloc[0]["purity"]
        else:
            purity = "not_available"

        # Add the purity to all SNVs as its own column.
        data["purity"] = purity
        data["path_to_purity_file"] = path_to_pcawg_purity_file

        data.to_csv(path_to_vcf, index=None, sep="\t")

    return path_to_vcf


def _parse_master_file_for_purity(purity_file: str):
    # Get the median after ESD, SD after ESD, and median before ESD.
    median_after_esd = None
    sd_after_esd = None
    median_before_esd = None

    with open(purity_file, "r") as f:
        purity_file_lines = f.readlines()
        for idx, row in enumerate(purity_file_lines):
            if "Median purity before ESD" in row:
                median_before_esd = float(
                    str(purity_file_lines[idx + 1]).strip())
            elif "Median purity after ESD" in row:
                median_after_esd = float(
                    str(purity_file_lines[idx + 1]).strip())
            elif "Standard deviation after ESD" in row:
                sd_after_esd = float(str(purity_file_lines[idx + 1]).strip())

    return median_after_esd, sd_after_esd, median_before_esd


@execution_time
def _add_strand_information(path: str, gencode_path: str):

    # Load in the data.
    data = pd.read_csv(path, delimiter="\t")
    if "strand" in list(data.columns):
        return path

    gencode_data = pd.read_csv(gencode_path, delimiter="\t")

    # Turn #chr in gencode_data into the same format as data.
    gencode_data["#chr"] = gencode_data["#chr"].str.replace(
        'chr', '').astype(str)

    # Create a dictionary off of the gencode data.
    gencode_dict = {}
    for _, row in gencode_data.iterrows():
        chrom = str(row["#chr"])
        gene = row["gene_name"]
        strand = row["strand"]

        if chrom not in gencode_dict:
            gencode_dict[chrom] = {}

        gencode_dict[chrom][gene] = strand

    # Create a list of the strands that should be appended to the data dataframe.
    strands = []
    for _, row in data.iterrows():
        chrom = str(row["#CHROM"])
        gene = row["GENE"]
        try:
            strands.append(str(gencode_dict[chrom][gene]))
        except:
            strands.append("unknown")

    # Define the new column "strand" to be the strands.
    data["strand"] = strands

    # Write the csv to the same path.
    data.to_csv(path, sep="\t", index=False)

    return path


@execution_time
def _add_ge_data_to_genes_and_transcription_factors(path_to_vcf_file, config):
    
    # Define the pid, which is generated from the .vcf file.
    pid = _get_pid_from_structured_vcf_path(
        path_to_vcf=path_to_vcf_file,
        only_pid=True
    )

    # Define variables from the configuration dictionary.
    dataset = config["pipeline"]["dataset"]
    path_to_rna_seq = config["rnaseq_data"][dataset]["path_to_rnaseq_dataframe"]
    path_to_metadata = config["pipeline"]["path_to_metadata"]
    gene_expression_measure = config["rnaseq_data"]["rnaseq_measurement"]
    name_of_tfbs_prediction_column = config["preprocessing_details"][
        "transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"]
 
    # Load in the dataframes.
    rna_seq = pd.read_csv(path_to_rna_seq, delimiter="\t",
                          index_col="Unnamed: 0")
    metadata = pd.read_csv(path_to_metadata)
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")

    # Check if the gene_expression_measure and the gene_expression_measure Z-Score is already in the dataframe.
    # If so, do nothing an return the path to the .vcf file. Else, continue with the rest of the code.
    if (gene_expression_measure in list(data.columns)) and (f"{gene_expression_measure}_Z_score" in list(data.columns)):
        print(f"{path_to_vcf_file} already has expression measurements!")
        return path_to_vcf_file

    # Create a dictionary such that the keys are the pids and the values are the cohort.
    # Ex: pid_to_cohort_dict["pid1"] = "cohort1"
    pid_to_cohort_dict = metadata.groupby('pid')['cohort'].first().to_dict()

    # Get the cohort of the patient.
    cohort_of_sample = pid_to_cohort_dict[pid]

    # Get the pids of the same cohort.
    pids_within_patient_cohort = list(metadata[(metadata["cancer_type"] == "tumor") & (
        metadata["cohort"] == cohort_of_sample)]["pid"])
    print(
        f"Number of pids within patient cohort: {len(pids_within_patient_cohort)}")

    # Remove those pids from pids_within_patient_cohort that are not within the rna_seq dataframe.
    pids_within_patient_cohort_with_rna_seq = list(
        set([pid for pid in pids_within_patient_cohort if pid in list(rna_seq.index)]))
    print(
        f"Number of pids within patient cohort with RNA-Seq: {len(pids_within_patient_cohort_with_rna_seq)}")

    #### (1) Add gene expression ####

    # Check if the pid has corresponding RNA-seq data available.
    if pid not in pids_within_patient_cohort_with_rna_seq:
        data[gene_expression_measure] = "not_available"
        data[f"{gene_expression_measure}_Z_score"] = "not_available"
    else:
        # Choose those pids of the rna_seq dataframe that are within pids_within_patient_cohort and subset the rna_seq dataframe.
        rna_seq_of_cohort_raw = rna_seq.loc[pids_within_patient_cohort_with_rna_seq]
        rna_seq_of_cohort_zscore = zscore(
            rna_seq_of_cohort_raw, axis=0).fillna(0)
        genes_with_rnaseq_data = list(rna_seq_of_cohort_raw.columns)

        # Get the RNA-seq column for the raw and zscore dataframes.
        if type(rna_seq_of_cohort_raw.loc[pid]) == pd.DataFrame:
            rna_seq_raw_col = pd.DataFrame(
                rna_seq_of_cohort_raw.loc[pid].iloc[0])
            rna_seq_zscore_col = pd.DataFrame(
                rna_seq_of_cohort_zscore.loc[pid].iloc[0])
        else:
            rna_seq_raw_col = pd.DataFrame(rna_seq_of_cohort_raw.loc[pid])
            rna_seq_zscore_col = pd.DataFrame(
                rna_seq_of_cohort_zscore.loc[pid])

        rna_seq_raw_col.columns = [f"{gene_expression_measure}"]
        rna_seq_zscore_col.columns = [f"{gene_expression_measure}_Z_score"]

        # Merge based on the gene names.
        data = pd.merge(data, rna_seq_raw_col, how='left',
                        left_on="GENE", right_index=True)
        data = pd.merge(data, rna_seq_zscore_col, how='left',
                        left_on="GENE", right_index=True)

    #### (2) Add transcription factor expression ####
    for idx, row in tqdm(data.iterrows(), total=data.shape[0]):
        transcription_factor_column = row[name_of_tfbs_prediction_column]

        new_transcription_factor_column = ""

        if transcription_factor_column == ".":
            new_transcription_factor_column = "."
        else:
            for entry in transcription_factor_column.split(";"):
                tf_name, binding_affinity, ref_seq, alt_seq = entry.split(",")
                try:
                    tf_raw = rna_seq_raw_col.loc[tf_name,
                                                 gene_expression_measure]
                    tf_zscore = rna_seq_zscore_col.loc[tf_name,
                                                       f"{gene_expression_measure}_Z_score"]
                    tf_log = np.log(tf_raw)
                except:
                    tf_raw = "not_available"
                    tf_zscore = "not_available"
                    tf_log = "not_available"
                new_transcription_factor_column += f"{tf_name},{binding_affinity},{ref_seq},{alt_seq},{tf_raw},{tf_zscore},{tf_log};"
            new_transcription_factor_column = new_transcription_factor_column[:-1]
        data.loc[idx,
                 f"{name_of_tfbs_prediction_column}(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"] = new_transcription_factor_column

    data.to_csv(
        path_to_vcf_file,
        sep="\t", index=False,
    )

    return path_to_vcf_file


def _get_sequence_context(
    chromosome: str,
    position: int,
    reference_path: str,
    context_length: int = 10,
    with_comma: bool = True,
    with_bracket: bool = False,
):
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"

    if not type(position) == int:
        position = int(position)

    # Load the reference genome
    ref = Fasta(reference_path)

    # Get the sequence using the chromosome and position
    sequence_context = ref[chromosome][position -
                                       context_length-1:position+context_length].seq

    # At the SNV location, replace the reference nucleotide with a comma if with_comma=True.
    if with_comma:
        sequence_context = list(sequence_context)
        sequence_context[context_length] = ","
        sequence_context = "".join(sequence_context)

    # At the SNV location, put brackets around the nucleotide if with_bracket=True.
    # This is used mainly for testing purposes.
    if with_bracket:
        sequence_context = list(sequence_context)
        sequence_context[context_length] = f"[{sequence_context[context_length]}]"
        sequence_context = "".join(sequence_context)

    return sequence_context.upper()


@execution_time
def _pcawg_add_sequence_context(
    path_to_vcf_file: str,
    path_to_hg19_fa_file: str,
):
    data = pd.read_csv(path_to_vcf_file, delimiter="\t")
    if "SEQUENCE_CONTEXT" not in list(data.columns):
        for idx, row in tqdm(data.iterrows(), total=data.shape[0]):
            chromosome = row["#CHROM"]
            position = row["POS"]
            sequence_context = _get_sequence_context(
                chromosome=str(chromosome),
                position=str(position),
                reference_path=path_to_hg19_fa_file
            )

            data.loc[idx, "SEQUENCE_CONTEXT"] = sequence_context

        data.to_csv(
            path_to_vcf_file,
            sep="\t", index=False,
        )

    return path_to_vcf_file

@execution_time
def _add_expression_to_gene_and_tf_of_prospective(
    path_to_sample_vcf: str,
    sample_config: dict,
    path_to_sample_metadata: str  = "/omics/groups/OE0436/internal/nabad/_final_results_29March2024/single_file/2023_week02_to_week11_metadata_with_updated_cohort.csv",
    path_to_prior_metadata: str  = "/omics/groups/OE0436/internal/nabad/_final_results_29March2024/MASTER_2022-10-16_13h03/master_metadata_2022-10-16_13h03.csv",
    path_to_prior_expression_data: str  = "/omics/groups/OE0436/internal/nabad/_final_results_29March2024/MASTER_2022-10-16_13h03/ge_data/raw_fpkm_dataframe.tsv",
): 
    expression_measurement = sample_config["rnaseq_data"]["rnaseq_measurement"]

    # Get the sample's dataframe.
    sample_dataframe = pd.read_csv(path_to_sample_vcf, delimiter="\t")

    # Ensure that the correct columns are in place within the sample's dataframe.
    name_of_gene_column = "GENE"
    name_of_jaspar_column = sample_config["preprocessing_details"]["transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"]
    assert name_of_gene_column in list(sample_dataframe.columns), f"{name_of_gene_column} is not in the columns of the dataframe."
    assert name_of_jaspar_column in list(sample_dataframe.columns), f"{name_of_jaspar_column} is not in the columns of the dataframe."

    # Get the cohort of the sample.
    sample_metadata = pd.read_csv(path_to_sample_metadata)
    pid = _get_pid_from_structured_vcf_path(path_to_sample_vcf, True)

    pid_metadata_subset = sample_metadata[sample_metadata["pid"] == pid]
    if pid_metadata_subset.shape[0] > 0:
        sample_cohort = pid_metadata_subset.iloc[0]["cohort"]
    else:
        sample_cohort = "Anderes"

    # Get the path to the sample's RNA-Seq file.
    if pid_metadata_subset.shape[0] > 0:
        sample_path_to_tsv = pid_metadata_subset.iloc[0]["path_to_tsv"]
        if not os.path.exists(sample_path_to_tsv):
            sample_path_to_tsv = "not_available"
    else:
        sample_path_to_tsv = "not_available"

    # Load in the prior metadata of the dataset to compare the zscore to.
    prior_metadata = pd.read_csv(path_to_prior_metadata)
    prior_pids_in_same_cohort_as_sample = list(prior_metadata[prior_metadata["cohort"] == sample_cohort]["pid"].unique())
    print(f"There are {len(prior_pids_in_same_cohort_as_sample)} retrospective pids in the same cohort.")

    # Load in the prior raw expression file.
    print(f"Loading in expression dataframe...")
    prior_raw_fpkm = pd.read_csv(path_to_prior_expression_data, delimiter="\t")
    prior_raw_fpkm = prior_raw_fpkm.rename({"Unnamed: 0": "pid"}, axis=1).set_index("pid")
    print(f"... done.")

    # Subset the prior raw expression dataframe to only include those within the same cohort as the sample.
    prior_pids_in_same_cohort_as_sample_in_index = [pid for pid in prior_pids_in_same_cohort_as_sample if pid in list(prior_raw_fpkm.index)]
    print(f"Of the {len(prior_pids_in_same_cohort_as_sample)} retrospective pids, {len(prior_pids_in_same_cohort_as_sample_in_index)} have expression data available within the expression dataframe.")

    cohort_prior_raw_fpkm = prior_raw_fpkm.loc[prior_pids_in_same_cohort_as_sample_in_index]

    # Load in the sample's RNA-Seq file.
    if os.path.exists(sample_path_to_tsv):
        samples_rnaseq = pd.read_csv(sample_path_to_tsv, delimiter="\t")

        # Create a dictionary from the sample's gene to it's expression_measurement.
        sample_expression_dictionary = {}
        for idx, row in samples_rnaseq.iterrows():
            gene_name = row["name"]
            sample_expression_dictionary[gene_name] = float(row[expression_measurement])
        print(f"There are total of {len(sample_expression_dictionary.keys())} genes with expression of the sample.")

        # Create a brand new row in the cohort_prior_raw_fpkm that adds the sample's expression.
        for gene_name in tqdm(cohort_prior_raw_fpkm.columns, desc="Adding the sample's expression:"):
            try:
                cohort_prior_raw_fpkm.loc[pid, gene_name] = sample_expression_dictionary[gene_name]
            except:
                cohort_prior_raw_fpkm.loc[pid, gene_name] = 0

        # Compute the z-score per gene.
        cohort_prior_zscore_fpkm = zscore(cohort_prior_raw_fpkm, axis=0).fillna(0)

        # Create a dictionary for the sample for gene/TF to expression.
        samples_expression_dictionary = {}

        for gene in tqdm(cohort_prior_zscore_fpkm.columns, desc="Creating raw/zscore expression dict for sample:"):
            samples_expression_dictionary[gene] = {
                "raw": cohort_prior_raw_fpkm.loc[pid, gene],
                "zscore": cohort_prior_zscore_fpkm.loc[pid, gene]
            }

    # At the gene level, add the raw expression and zscore.
    for idx, row in tqdm(sample_dataframe.iterrows(), desc = "1. Adding gene-level expression", total=sample_dataframe.shape[0]):
        gene = row["GENE"]

        if os.path.exists(sample_path_to_tsv):
            sample_dataframe.loc[idx, expression_measurement] = samples_expression_dictionary[gene]["raw"] if gene_name in samples_expression_dictionary else "not_available"
            sample_dataframe.loc[idx, f"{expression_measurement}_Z_score"] = samples_expression_dictionary[gene]["zscore"]  if gene_name in samples_expression_dictionary else "not_available"
        else:
            sample_dataframe.loc[idx, expression_measurement] = "not_available"
            sample_dataframe.loc[idx, f"{expression_measurement}_Z_score"] = "not_available"

    # At the TF level, add the raw expression, zscore and log.
    for idx, row in tqdm(sample_dataframe.iterrows(), desc = "2. Add TF-level expression", total=sample_dataframe.shape[0]):
        transcription_factor_column = row[name_of_jaspar_column]

        new_transcription_factor_column = ""

        if transcription_factor_column == ".":
            new_transcription_factor_column = "."
        else:
            for entry in transcription_factor_column.split(";"):
                tf_name, binding_affinity, ref_seq, alt_seq = entry.split(",")
                if os.path.exists(sample_path_to_tsv):
                    if tf_name in samples_expression_dictionary:
                        tf_raw = str(float(samples_expression_dictionary[tf_name]["raw"])) if samples_expression_dictionary[tf_name]["raw"] != "not_available" else "not_available"
                        tf_zscore = str(float(samples_expression_dictionary[tf_name]["zscore"])) if samples_expression_dictionary[tf_name]["zscore"] != "not_available" else "not_available"
                        tf_log = str(np.log(float(samples_expression_dictionary[tf_name]["raw"]))) if samples_expression_dictionary[tf_name]["raw"] != "not_available" else "not_available"
                    else:
                        tf_raw, tf_zscore, tf_log = "not_available", "not_available", "not_available"
                else:
                    tf_raw, tf_zscore, tf_log = "not_available", "not_available", "not_available"

                new_transcription_factor_column += f"{tf_name},{binding_affinity},{ref_seq},{alt_seq},{tf_raw},{tf_zscore},{tf_log};"

            new_transcription_factor_column = new_transcription_factor_column[:-1]
        sample_dataframe.loc[idx, f"{name_of_jaspar_column}(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"] = new_transcription_factor_column
        
    sample_dataframe.to_csv(
        path_to_sample_vcf,
        sep="\t", index=False,
    )


@execution_time
def main(
    path_to_vcf_file: str,
    config: dict
):
    dataset = config["pipeline"]["dataset"]
    prospective = config["pipeline"]["prospective"]

    # Rename the .vcf file to have the suffix.
    suffix_to_append_to_vcf = config["preprocessing_details"]["suffix_to_append_to_vcf"]
    preprocessed_filename = path_to_vcf_file.replace(
        ".vcf", suffix_to_append_to_vcf)

    # Copy the contents of the .vcf file to the new file location.
    shutil.copy(
        path_to_vcf_file,
        preprocessed_filename
    )

    if dataset == "pcawg":
        preprocessed_filename = _pcawg_rename_hugo_symbol_to_gene(
            preprocessed_filename)

        path_to_hg19_reference = config["preprocessing_details"]["pcawg"]["path_to_genome_reference_fa_file"]
        preprocessed_filename = _pcawg_add_sequence_context(
            preprocessed_filename, path_to_hg19_reference)

    preprocessed_filename = _fix_problematic_genes(preprocessed_filename)

    preprocessed_filename = _remove_lines_with_no_gene(preprocessed_filename)

    preprocessed_filename = _remove_indels(preprocessed_filename)

    preprocessed_filename = _remove_snp_and_germline(preprocessed_filename)

    preprocessed_filename = _add_purity_information(
        preprocessed_filename, config)

    preprocessed_filename = _add_vaf_information_to_final_vcf_files(
        preprocessed_filename, config["pipeline"]["dataset"])

    preprocessed_filename = _add_strand_information(
        preprocessed_filename, config["additional_files"]["tss_reference_file"])

    preprocessed_filename = run_fimo(
        path_to_vcf_file=preprocessed_filename,
        path_to_fimo=config["preprocessing_details"]["transcription_factor_prediction"]["path_to_fimo_executable"],
        path_to_jaspar_database=config["preprocessing_details"][
            "transcription_factor_prediction"]["path_to_jaspar_database"],
        name_of_jaspar_column=config["preprocessing_details"][
            "transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"],
        save_tfbs_dict=config["preprocessing_details"]["transcription_factor_prediction"],
    )

    preprocessed_filename = _fix_transcription_factor_names(
        preprocessed_filename,
        config["preprocessing_details"]["transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"]
    )
    if prospective:
        preprocessed_filename = _add_expression_to_gene_and_tf_of_prospective(
            path_to_sample_vcf = preprocessed_filename,
            sample_config = config,
            path_to_sample_metadata = config["pipeline"]["path_to_metadata"],
            path_to_prior_metadata  = "/omics/groups/OE0436/internal/nabad/_final_results_29March2024/MASTER_2022-10-16_13h03/master_metadata_2022-10-16_13h03.csv",
            path_to_prior_expression_data = "/omics/groups/OE0436/internal/nabad/_final_results_29March2024/MASTER_2022-10-16_13h03/ge_data/raw_fpkm_dataframe.tsv",
        )
    else:
        preprocessed_filename = _add_ge_data_to_genes_and_transcription_factors(
            preprocessed_filename,
            config
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
        name_of_pipeline_step = "preprocessing",
        start_time = start_time,
        end_time = end_time,
        name_of_timing_tracker_file = "timing_tracker.json",
    )
    
    print("~~Finished with preprocessing~~")
