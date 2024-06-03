import os
from zipfile import ZipFile
from shutil import rmtree
from pipeline.args import GeneralArgs
from pipeline.propagation_routines import perform_propagation
from pipeline.pathway_enrichment import perform_enrichment
from pipeline.utils import read_temp_scores, process_condition
from pipeline.visualization_tools import print_aggregated_pathway_information, plot_pathways_mean_scores

def run_pipeline(run_propagation=True, run_enrichment=True):
    general_args = GeneralArgs(run_propagation=run_propagation, alpha=0.1)

    test_file_paths = [os.path.join(general_args.input_dir, file) for file in os.listdir(general_args.input_dir) if file.endswith('.xlsx')]
    test_name_list = [os.path.splitext(file)[0] for file in os.listdir(general_args.input_dir) if file.endswith('.xlsx')]

    for test_name in test_name_list:
        if run_propagation:
            print(f"Running propagation on {test_name}")
            perform_propagation(test_name, general_args)

        if run_enrichment:
            print(f"Running enrichment on {test_name}")
            perform_enrichment(test_name, general_args)
            print("-----------------------------------------")

    print("Finished enrichment")

    condition_files = [os.path.join(general_args.temp_output_folder, file) for file in os.listdir(general_args.temp_output_folder)]
    all_pathways = {}

    for condition_file in condition_files:
        enriched_pathway_dict = read_temp_scores(condition_file)
        for pathway in enriched_pathway_dict:
            if pathway not in all_pathways:
                all_pathways[pathway] = {}

    for condition_file, experiment_file in zip(condition_files, test_file_paths):
        process_condition(condition_file, experiment_file, general_args.pathway_file_dir, all_pathways)

    aggregated_text_file = print_aggregated_pathway_information(general_args, all_pathways)
    plot_file = plot_pathways_mean_scores(general_args, all_pathways)

    if os.path.exists(general_args.temp_output_folder):
        rmtree(general_args.temp_output_folder)

    result_zip_path = os.path.join(general_args.output_dir, "results.zip")
    with ZipFile(result_zip_path, 'w') as result_zip:
        result_zip.write(aggregated_text_file, arcname=os.path.basename(aggregated_text_file))
        result_zip.write(plot_file, arcname=os.path.basename(plot_file))

    return result_zip_path
