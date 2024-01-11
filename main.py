import time
from args import NewEnrichTask, NewGeneralArgs
from statistic_methods import kolmogorov_smirnov_test
from pathway_enrichment import run

def perform_enrichment():
    # run enrichment
    general_args = NewGeneralArgs()
    print("running enrichment")
    task1 = NewEnrichTask(experiment_name='roded_T_v_N', statistic_test=kolmogorov_smirnov_test)
    print('running')
    run(task1, general_args)


def main(run_rop, run_enrichment):
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
    if run_enrichment:
        perform_enrichment()


if __name__ == '__main__':
    start = time.time()

    run_propagation_flag = False
    run_enrichment_flag = True

    main(run_propagation_flag, run_enrichment_flag)

    end = time.time()

    print("Time elapsed: {} seconds".format(end - start))