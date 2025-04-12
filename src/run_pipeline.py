import warnings
import json
import time
import argparse

warnings.simplefilter(action='ignore')

from pipeline.preprocessing import run_preprocessing_on_all_paths
from pipeline.filters.ge_filter import run_ge_on_all_paths
from pipeline.filters.promoter_filter import run_promoter_on_all_paths
from pipeline.filters.motif_and_tf_filter import run_all_motif_and_tf_expression_on_all_paths
from pipeline.annotations.cgc_list import add_cgc_to_all_paths
from pipeline.annotations.chromhmm import add_chromhmm_to_all_paths
from pipeline.annotations.recurrence import add_recurrence_to_all_paths
from pipeline.postprocessing import run_postprocessing_on_all_paths

def run_initial_steps(config_path):
    """
    Runs the initial pipeline steps based on the provided configuration.

    Parameters
    ----------
    config : str
        The pipeline configuration loaded from a JSON file.

    Raises
    ------
    ValueError
        If an invalid initial step is encountered.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
    
    for step in config["pipeline"]["initial_steps"]:
        print(f"####### {step.upper()} #######")
        if step.lower() == "initial_structure":
            pass  # Placeholder for future implementation
        elif step.lower() == "preprocessing":
            run_preprocessing_on_all_paths.main(config_path)
        else:
            raise ValueError(f"Invalid initial step: {step}")

def run_filters(config_path):
    """
    Runs filtering steps in the pipeline.

    Parameters
    ----------
    config : str
        The pipeline configuration loaded from a JSON file.

    Raises
    ------
    ValueError
        If an invalid filter step is encountered.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
    
    for filter_step in config["pipeline"]["filter_order"]:
        print(f"####### {filter_step.upper()} #######")
        if filter_step.lower() == "promoter_filter":
            run_promoter_on_all_paths.main(config_path)
        elif filter_step.lower() == "ge_filter":
            run_ge_on_all_paths.main(config_path)
        elif filter_step.lower() == "motif_and_tf_expression_filter":
            run_all_motif_and_tf_expression_on_all_paths.main(config_path)
        else:
            raise ValueError(f"Invalid filter: {filter_step}")

def run_annotations(config_path):
    """
    Runs annotation steps in the pipeline.

    Parameters
    ----------
    config : str
        The pipeline configuration loaded from a JSON file.

    Raises
    ------
    ValueError
        If an invalid annotation step is encountered.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
    
    for annotation in config["pipeline"]["additional_annotations"]:
        print(f"####### {annotation.upper()} #######")
        if annotation.lower() == "cgc_annotation":
            add_cgc_to_all_paths.main(config_path)
        elif annotation.lower() == "chromhmm_annotation":
            add_chromhmm_to_all_paths.main(config_path)
        elif annotation.lower() == "recurrence_annotation":
            add_recurrence_to_all_paths.main(config_path)
        else:
            raise ValueError(f"Invalid annotation: {annotation}")

def run_post_pipeline_steps(config_path):
    """
    Runs post-processing steps in the pipeline.

    Parameters
    ----------
    config : str
        The pipeline configuration loaded from a JSON file.

    Raises
    ------
    ValueError
        If an invalid post-processing step is encountered.
    """
    with open(config_path, "r") as f:
        config = json.load(f)
        
    for step in config["pipeline"]["post_pipeline_steps"]:
        print(f"####### {step.upper()} #######")
        if step.lower() == "postprocessing":
            run_postprocessing_on_all_paths.main(config_path)
        else:
            raise ValueError(f"Invalid post-pipeline step: {step}")

def main(config_path):
    """
    Runs the full pipeline based on the provided configuration file.

    Parameters
    ----------
    config_path : str
        Path to the configuration JSON file.
    """
    start_time = time.time()
    
    print("\nStarting Pipeline Execution...\n")

    run_initial_steps(config_path)
    run_filters(config_path)
    run_annotations(config_path)
    run_post_pipeline_steps(config_path)
    
    print("\nPipeline execution complete. Total time: {:.2f} seconds".format(time.time() - start_time))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the REMIND-Cancer Pipeline")
    parser.add_argument("-c", "--config", required=True, help="Path to the configuration file.")
    
    args = parser.parse_args()
    main(args.config)
