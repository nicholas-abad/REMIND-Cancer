import warnings
warnings.simplefilter(action='ignore')

import json
import time
from optparse import OptionParser

# from pipeline_setup.folder_structure import create_initial_structure

from pipeline.preprocessing import run_preprocessing_on_all_paths
from pipeline.filters.ge_filter import run_ge_on_all_paths
from pipeline.filters.promoter_filter import run_promoter_on_all_paths
from pipeline.filters.motif_and_tf_filter import run_all_motif_and_tf_expression_on_all_paths
from pipeline.annotations.cgc_list import add_cgc_to_all_paths
from pipeline.annotations.chromhmm import add_chromhmm_to_all_paths
from pipeline.annotations.recurrence import add_recurrence_to_all_paths
from pipeline.postprocessing import run_postprocessing_on_all_paths


if __name__ == "__main__":
    starting_time = time.time()
    parser = OptionParser()
    parser.add_option(
        "-c",
        "--config",
        action="store",
        type="str",
        dest="config",
        help="Path to the configuration file.\n",
    )

    (options, args) = parser.parse_args()
    
    # Turn configuration file into a dictionary.
    with open(options.config, "r") as f:
        config = json.load(f)

    # Run intial steps.
    initial_steps = config["pipeline"]["initial_steps"]
    for step in initial_steps:
        if step.lower() == "initial_structure":
            print("####### initial structure #######")
            # create_initial_structure.main(options.config)
        elif step.lower() == "preprocessing":
            print("####### preprocessing #######")
            # run_preprocessing_on_all_paths.main(options.config)
        else:
            assert False, f"The step '{step}' does not exist."

    # Run pipeline.
    filter_order = config["pipeline"]["filter_order"]
    for filter in filter_order:
        if filter.lower() == "promoter_filter":
            print("####### promoter_filter #######")
            # run_promoter_on_all_paths.main(options.config)
        elif filter.lower() == "ge_filter":
            print(" ####### GE FILTER ########")
            run_ge_on_all_paths.main(options.config)
        elif filter.lower() == "motif_and_tf_expression_filter":
            print(" ####### MOTIF AND TF EXPRESSION FILTER ########")
            run_all_motif_and_tf_expression_on_all_paths.main(options.config)
        else:
            assert False, f"The filter '{filter}' does not exist."
    
    # Add additional annotations.
    annotations_order = config["pipeline"]["additional_annotations"]
    for annotation in annotations_order:
        if annotation.lower() == "cgc_annotation":
            print(" ####### CGC ANNOTATION ########")
            add_cgc_to_all_paths.main(options.config)
            
        elif annotation.lower() == "chromhmm_annotation":
            print(" ####### Open Chromatin / ChromHMM ANNOTATION ########")
            add_chromhmm_to_all_paths.main(options.config)
            
        elif annotation.lower() == "recurrence_annotation":
            print(" ####### RECURRENCE ANNOTATION ########")
            add_recurrence_to_all_paths.main(options.config)
        else:
            assert False, f"The annotation '{annotation}' does not exist."

    # Run post-pipeline steps.
    post_pipeline_steps = config["pipeline"]["post_pipeline_steps"]
    for step in post_pipeline_steps:
        print(f"Step: {step}")
        if step.lower() == "postprocessing":
            print(" ####### POST PROCESSING ########")
            run_postprocessing_on_all_paths.main(options.config)
        else:
            assert False, f"The step '{step}' does not exist."

    print("Complete pipeline time: ", time.time() - starting_time)
