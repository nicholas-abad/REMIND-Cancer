import json
import pandas as pd
from tqdm import tqdm
from optparse import OptionParser
import os

# For proper import structure:
import sys
path_to_remind_cancer_folder = os.getcwd()
if path_to_remind_cancer_folder not in sys.path:
    sys.path.append(path_to_remind_cancer_folder)

from src.pipeline.general_helper_functions import execution_time

def _get_tfs_passing_motif_and_tf_expression_filters(
    transcription_factor_entries: str,
    tfbs_destruction_threshold: float=0.08,
    tfbs_creation_threshold: float=11,
    tf_expression_measure_to_filter: str="raw",
    tf_expression_threshold: float=0  
)-> str:
    """Given a TF entry, keep the _entire_ entry if one or more transcription factors pass both the 
    motif filter _and_ the transcription factor expression filter. A transcription factor passes the
    motif filter if its binding affinity is less than or equal to the destruction threshold _or_ if
    its greater than or equal to the creation threshold. A transcription factor pass the TF expression
    filter if the tf_expression_measure_to_filter (i.e. raw) is strictly greater than the tf_expression_threshold.

    If there is not at least one TF within all of the entries or if transcription_factor_entries is ".",
    the returned string is ".".

    Example 1:
        Input:
            "some_tf,11,CCGGGTCCCCGG,.,10,10,30;some_tf,11,CCGGGTCCCCGG,.,not_available,not_available,not_available"
        Output:
            "some_tf,11,CCGGGTCCCCGG,.,10,10,30;some_tf,11,CCGGGTCCCCGG,.,not_available,not_available,not_available"

    Example 2:
        Input:
            "some_tf,11,CCGGGTCCCCGG,.,0,10,30;some_tf,11,CCGGGTCCCCGG,.,not_available,not_available,not_available"
        Output:
            "."


    Parameters
    ----------
    transcription_factor_entries : str
        Transcription factor entries of a single cohort in the format "tf_name,binding_affinity,seq1,seq2,raw,zscore,log".
        Ex: "TEST_TF1,11,CCGGGTCCCCGG,.,10,10,30;TEST_TF2,0,CCGGGTCCCCGG,.,0,-1,30"
    tfbs_destruction_threshold : float, optional
        Threshold for defining a destroyed transcription factor binding site, by default 0.01.
    tfbs_creation_threshold : float, optional
        Threshold for defining a created transcription factor binding site, by default 11.
    tf_expression_measure_to_filter : str, optional
       Type of score to compare for the TF Expression Filter, by default "raw".
       This can be one of ["raw", "zscore", "z-score", or "log"]
    tf_expression_threshold : float, optional
        Threshold for defining the activity of a transcription factor, by default 0.

    Returns
    -------
    str
        Transcription factor entries that pass both the Motif Filter _and_ the TF Expression Filter.
    """
    assert tf_expression_measure_to_filter in ["raw", "zscore", "z-score", "log"]
    
    if transcription_factor_entries == ".":
        return "."
    
    at_least_one_tf_passes_motif_and_tf_exp_filters = False
    
    # Iterate through each transcription factor in transcription_factor_entries.
    for entry in transcription_factor_entries.split(";"):
        _, binding_affinity, _, _, raw, zscore, log = entry.split(",")
        
        if tf_expression_measure_to_filter == "raw":
            expression = raw
        elif tf_expression_measure_to_filter == "log":
            expression = log
        else:
            expression = zscore
            
        if expression == "not_available":
            continue
            
        # Check if this specific TF passes the motif and TF Expression filter.
        passes_motif_filter = (float(binding_affinity) <= tfbs_destruction_threshold) or (float(binding_affinity) >= tfbs_creation_threshold)
        passes_tf_expression_filter = float(expression) > tf_expression_threshold
        
        if passes_motif_filter and passes_tf_expression_filter:
            at_least_one_tf_passes_motif_and_tf_exp_filters = True
        
    if at_least_one_tf_passes_motif_and_tf_exp_filters:
        return transcription_factor_entries
    else:
        return "."

def _add_number_of_created_and_destroyed_tfs_passing_motif_and_tf_exp_filter(
    dataframe: pd.DataFrame,
    name_of_jaspar_column_with_ge_data: str="JASPAR2020_CORE_vertebrates_non_redundant(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)",
    tfbs_creation_threshold: float=11,
    tfbs_destruction_threshold: float=0.09,
    tf_expression_measure_to_filter: str="raw",
    tf_expression_threshold: float=0,
):
    assert tf_expression_measure_to_filter in ["raw", "zscore", "log"], "tf_expression_measure_to_filter must be 'raw', 'zscore', or 'log'."
    
    
    for idx, row in tqdm(dataframe.iterrows(), total=dataframe.shape[0], desc="Adding created/destroyed TF information"):
        tfs = row[name_of_jaspar_column_with_ge_data]
        
        tfs_created_and_passing_tf_expression_threshold = []
        tfs_destroyed_and_passing_tf_expression_threshold = []
        remaining_tfs = []
        print(f"TFs: {tfs}")
        for tf in tfs.split(";"):
            tf_name, binding_affinity, seq1, seq2, raw, zscore, log = tf.split(",")
            if raw == "not_available":
                continue
                
            if tf_expression_measure_to_filter == "raw":
                expression_to_use = float(raw)
            elif tf_expression_measure_to_filter == "zscore":
                expression_to_use = float(zscore)
            else:
                expression_to_use = float(log)
                
            binding_affinity = float(binding_affinity)
                        
            if (binding_affinity <= tfbs_destruction_threshold) and (expression_to_use > tf_expression_threshold):
                tfs_destroyed_and_passing_tf_expression_threshold.append(tf_name)
            elif (binding_affinity >= tfbs_creation_threshold) and (expression_to_use > tf_expression_threshold):
                tfs_created_and_passing_tf_expression_threshold.append(tf_name)
            else:
                remaining_tfs.append(tf_name)

        assert len(tfs_created_and_passing_tf_expression_threshold) + len(tfs_destroyed_and_passing_tf_expression_threshold) > 0

        dataframe.loc[idx, "created_tfs_passing_tf_expression_threshold"] = ",".join(tfs_created_and_passing_tf_expression_threshold) if len(tfs_created_and_passing_tf_expression_threshold) > 0 else "None"
        dataframe.loc[idx, "destroyed_tfs_passing_tf_expression_threshold"] = ",".join(tfs_destroyed_and_passing_tf_expression_threshold) if len(tfs_destroyed_and_passing_tf_expression_threshold) > 0 else "None"
        dataframe.loc[idx, "remaining_tfs"] = ",".join(remaining_tfs) if len(remaining_tfs) > 0 else "None"
        
        dataframe.loc[idx, "num_created_tfs_passing_tf_expression_threshold"] = len(tfs_created_and_passing_tf_expression_threshold)
        dataframe.loc[idx, "num_destroyed_tfs_passing_tf_expression_threshold"] = len(tfs_destroyed_and_passing_tf_expression_threshold)
        dataframe.loc[idx, "num_remaining_tfs"] = len(remaining_tfs)

    return dataframe

@execution_time
def _motif_and_tf_expression_filter(
    path_to_single_vcf: str,
    tfbs_destruction_threshold: float=0.08,
    tfbs_creation_threshold: float=11,
    tf_expression_measure_to_filter: str="raw",
    tf_expression_threshold: float=0,
    suffix_to_append_to_vcf: str="_after_combined_motif_and_tf_filter.vcf",
    name_of_jaspar_column_with_ge_data: str="JASPAR2020_CORE_vertebrates_non_redundant(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"
) -> None:
    """Run the motif and transcription factor expression filters.

    Parameters
    ----------
    path_to_single_vcf : str
        Path to a patient's .vcf.
    tfbs_destruction_threshold : float, optional
        Threshold for defining a destroyed transcription factor binding site, by default 0.01.
    tfbs_creation_threshold : float, optional
        Threshold for defining a created transcription factor binding site, by default 11.
    tf_expression_measure_to_filter : str, optional
       Type of score to compare for the TF Expression Filter, by default "raw".
       This can be one of ["raw", "zscore", "z-score", or "log"]
    tf_expression_threshold : float, optional
        Threshold for defining the activity of a transcription factor, by default 0.
    suffix_to_append_to_vcf : str, optional
        Tag to add to the end of the .vcf file, by default "_after_combined_motif_and_tf_filter.vcf".
    """
    # Load in the dataframe.
    data = pd.read_csv(path_to_single_vcf, delimiter="\t")
    
    # Get the proper JASPAR column name, which contains the predicted transcription factors.

    indices_to_keep = []
    
    for idx, row in tqdm(data.iterrows(), total=data.shape[0]):
        # If there is at least 1 TF in at least 1 cohort that passes both motif and tf_expression filter, keep. Else, delete the row.
        keep_row = False

        transcription_factor_entries = row[name_of_jaspar_column_with_ge_data]

        # Get the transcription factor entries where at least one TF passes both filters.
        new_transcription_factor_entries = _get_tfs_passing_motif_and_tf_expression_filters(
            transcription_factor_entries = transcription_factor_entries,
            tfbs_destruction_threshold = tfbs_destruction_threshold,
            tfbs_creation_threshold = tfbs_creation_threshold,
            tf_expression_measure_to_filter = tf_expression_measure_to_filter,
            tf_expression_threshold = tf_expression_threshold
        )

        # Replace the old TF entry with the new one.
        data.loc[idx, name_of_jaspar_column_with_ge_data] = new_transcription_factor_entries

        # Change keep_row to True if the new TF entries are not empty.
        if new_transcription_factor_entries != ".":
            keep_row = True
                
        if keep_row:
            indices_to_keep.append(idx)
            

    data = data.iloc[indices_to_keep]
    
    if data.shape[0] > 0:
        data = _add_number_of_created_and_destroyed_tfs_passing_motif_and_tf_exp_filter(
            dataframe = data,
            name_of_jaspar_column_with_ge_data = name_of_jaspar_column_with_ge_data,
            tfbs_creation_threshold = tfbs_creation_threshold,
            tfbs_destruction_threshold = tfbs_destruction_threshold,
            tf_expression_measure_to_filter = tf_expression_measure_to_filter,
            tf_expression_threshold = tf_expression_threshold,
        )
        
    
        data.to_csv(
            path_to_single_vcf.replace(".vcf", suffix_to_append_to_vcf),
            index=False, sep="\t"
        )

@execution_time
def main(path_to_vcf_file: str, config: dict):
    
    suffix_to_append_to_vcf = config["motif_and_tf_expression_filter"]["suffix_to_append_to_vcf"]
    tfbs_creation_threshold = config["motif_and_tf_expression_filter"]["tfbs_creation_threshold"]
    tfbs_destruction_threshold = config["motif_and_tf_expression_filter"]["tfbs_destruction_threshold"]
    tf_expression_measure_to_filter = config["motif_and_tf_expression_filter"]["tf_expression_measure_to_filter"]
    tf_expression_threshold = config["motif_and_tf_expression_filter"]["tf_expression_threshold"]
    
    # Get the name of the jaspar column with the 
    name_of_jaspar_column_with_ge_data = config["preprocessing_details"]["transcription_factor_prediction"]["name_of_tfbs_prediction_column_to_add"] + \
        "(tf_name,binding_affinity,seq1,seq2,raw,zscore,log)"
    
    _motif_and_tf_expression_filter(
        path_to_single_vcf = path_to_vcf_file,
        tfbs_destruction_threshold = tfbs_destruction_threshold,
        tfbs_creation_threshold = tfbs_creation_threshold,
        tf_expression_measure_to_filter = tf_expression_measure_to_filter,
        tf_expression_threshold = tf_expression_threshold,
        suffix_to_append_to_vcf = suffix_to_append_to_vcf,
        name_of_jaspar_column_with_ge_data = name_of_jaspar_column_with_ge_data
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
    
    print("~~ Finished with Motif and TF Expression filter ~~")
    