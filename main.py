import time
import shutil
from os import path, listdir
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import read_temp_scores, process_condition
from visualization_tools import print_aggregated_pathway_information, plot_pathways_mean_scores


def main(run_propagation=True, run_enrichment=True):
    """
    Main function to execute propagation and enrichment analysis based on specified flags.
    This function initializes tasks for propagation and enrichment and executes them based on the
    provided flags. It serves as the entry point for running the gene score propagation and enrichment analysis pipeline.
    Parameters:
    - run_propagation (bool): Flag to determine whether to run propagation (default: True).
    - run_enrichment (bool): Flag to determine whether to run enrichment analysis (default: True).
    Returns:
    - None: This function orchestrates the execution of other functions but does not return a value.
    """
    general_args = GeneralArgs(run_propagation=run_propagation)

    # Get a list of all .xlsx files in the input directory
    test_file_paths = [path.join(general_args.input_dir, file) for file in listdir(general_args.input_dir) if file.endswith('.xlsx')]

    # Get a list of just the names of the files without the extensions
    test_name_list = [path.splitext(file)[0] for file in listdir(general_args.input_dir) if file.endswith('.xlsx')]

    for test_name in test_name_list:
        if run_propagation:
            print(f"Running propagation on {test_name}")
            perform_propagation(test_name, general_args)

        if run_enrichment:
            print(f"Running enrichment on {test_name}")
            perform_enrichment(test_name, general_args)

    print("finished enrichment")

    # Get the list of condition files
    condition_files = [path.join(general_args.temp_output_folder, file) for file in listdir(general_args.temp_output_folder)]

    all_pathways = {}

    # Load enriched pathways from files into a dictionary for further processing
    for condition_file in condition_files:
        enriched_pathway_dict = read_temp_scores(condition_file)
        for pathway in enriched_pathway_dict.keys():
            if pathway not in all_pathways:
                all_pathways[pathway] = {}

    # Process conditions and aggregate data
    for condition_file, experiment_file in zip(condition_files, test_file_paths):
        process_condition(condition_file, experiment_file, general_args.pathway_file_dir, all_pathways)

    # Output aggregated pathway information to a text file
    print_aggregated_pathway_information(general_args.output_dir, general_args.Experiment_name, all_pathways)

    # Visualize mean scores of pathways across all conditions
    plot_pathways_mean_scores(general_args.output_dir, general_args.Experiment_name, all_pathways)

    # Clean up temporary output folder if it exists
    if path.exists(general_args.temp_output_folder):
        shutil.rmtree(general_args.temp_output_folder)

if __name__ == '__main__':
    start = time.time()
    # Set these flags to control the tasks to run
    run_propagation_flag = False
    run_enrichment_flag = True

    main(run_propagation=run_propagation_flag, run_enrichment=run_enrichment_flag)

    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))
