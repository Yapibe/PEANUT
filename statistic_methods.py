import numpy as np
from scipy.stats import ranksums, rankdata, ttest_ind, ks_2samp, ttest_rel,norm, hypergeom
from tqdm import tqdm
from math import sqrt, fabs


class StatResults:
    """
    Class to store the results of a statistical test.
    Attributes:
        p_value (float): The p-value of the test.
        directionality (bool): Indicates the direction of the difference (True if experiment > elements).
        name (str): Name of the statistical test.
    """

    def __init__(self, p_value=None, directionality=None, name=None):
        self.p_value = p_value
        self.directionality = directionality
        self.name = name


def wilcoxon_rank_sums_test(experiment_scores, control_scores, alternative='two-sided'):
    """
    Perform the Wilcoxon rank-sum test to compare two independent samples.

    Parameters:
    - experiment_scores (list): Scores from the experimental group.
    - control_scores (list): Scores from the control group.
    - alternative (str): Defines the alternative hypothesis ('two-sided', 'less', 'greater').

    Returns:
    float: The P-value from the Wilcoxon rank-sum test.
    """
    from scipy.stats import ranksums
    p_vals = ranksums(experiment_scores, control_scores, alternative=alternative).pvalue
    return p_vals


# Wilcoxon signed-rank test function
def wilcoxon_test(real_scores, shuffled_scores):
    # Wilcoxon signed-rank test logic
    differences = real_scores - shuffled_scores
    signed_ranks = rankdata(np.abs(differences)) * np.sign(differences)
    return np.sum(signed_ranks)


def sign_test(real_scores, shuffled_scores):
    differences = real_scores - shuffled_scores
    positive_count = np.sum(differences > 0)
    negative_count = np.sum(differences < 0)

    if positive_count > negative_count:
        return 1  # Indicate majority of differences are positive
    elif negative_count > positive_count:
        return -1  # Indicate majority of differences are negative
    else:
        return 0  # Equal number of positive and negative differences, or all zeros



def students_t_test(experiment_scores, elements_scores):
    """
    Performs Student's t-test or Welch's t-test to compare two sets of scores.

    Args:
        experiment_scores (array_like): Array of scores from the experiment group.
        elements_scores (array_like): Array of scores from the control group.

    Returns:
        StatResults: Object containing the p-value, directionality, and name of the test.
    """
    students_t_test.name = 'Students_t_test'

    # Calculate standard deviations
    std1 = np.std(experiment_scores, ddof=1)
    std2 = np.std(elements_scores, ddof=1)

    # Determine if variances are similar
    similar_variances = min(std1, std2) * 2 >= max(std1, std2)

    # Perform the appropriate t-test
    t_stat, p_vals = ttest_ind(experiment_scores, elements_scores, equal_var=similar_variances)
    direction = np.mean(experiment_scores) > np.mean(elements_scores)

    return StatResults(p_value=p_vals, directionality=direction, name=students_t_test.name)


def bh_correction(p_values):
    """
    Applies the Benjamini-Hochberg procedure for controlling the false discovery rate in multiple hypothesis testing.
    Args:
        p_values (array_like): Array of p-values obtained from multiple statistical tests.
    Returns:
        numpy.ndarray: Array of adjusted p-values.
    """
    p_vals_rank = rankdata(p_values, 'max') - 1
    p_vals_rank_ord = rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) / (p_vals_rank + 1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_vals_rank]
    return adj_p_vals


def calculate_empirical_p_values(real_data, pathways, test_function, num_simulations=1000):
    empirical_p_values = {}

    # Wrap the pathways iteration with tqdm for progress tracking
    for pathway_name, pathway_genes in tqdm(pathways.items(), desc='Processing Pathways'):
        if pathway_name == "WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV" or pathway_name == "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER":
            print("synaptic pathway")
        count_negative_t, count_positive_t = 0, 0
        pathway_mask = real_data['GeneID'].isin(pathway_genes)
        pathway_real_scores = real_data[pathway_mask]['Score']

        for i in range(1, num_simulations + 1):
            shuffled_score_column = f'Shuffled_Score_{i}'
            pathway_shuffled_scores = real_data[pathway_mask][shuffled_score_column]
            T = test_function(pathway_real_scores, pathway_shuffled_scores)
            count_negative_t += (T < 0)
            count_positive_t += (T > 0)

        min_count = min(count_negative_t, count_positive_t)
        empirical_p_values[pathway_name] = 2 * (min_count + 1) / (num_simulations + 1)

    return empirical_p_values


# Paired sample t-test function
def t_test(real_scores, shuffled_scores):
    t_statistic, _ = paired_sample_t_test(real_scores, shuffled_scores)
    return t_statistic


def paired_sample_t_test(real_scores, shuffled_scores):
    # Ensure that real_scores and shuffled_scores are NumPy arrays of type float
    real_scores = np.asarray(real_scores, dtype=float)
    shuffled_scores = np.asarray(shuffled_scores, dtype=float)

    # Perform the paired sample t-test
    t_statistic, p_value = ttest_rel(real_scores, shuffled_scores)
    return t_statistic


def kolmogorov_smirnov_test(experiment_scores, control_scores):
    """
    Perform the Kolmogorov-Smirnov test to compare two samples.

    Parameters:
    - experiment_scores (list): Scores from the experimental group.
    - control_scores (list): Scores from the control group.

    Returns:
    float: The P-value from the KS test indicating statistical difference.
    """
    # Convert lists to numpy arrays for efficient operations
    experiment_scores = np.sort(experiment_scores)
    control_scores = np.sort(control_scores)

    # Calculate the length of each sample
    en1 = len(experiment_scores)
    en2 = len(control_scores)

    # Combine the scores and compute cumulative distribution functions
    data_all = np.concatenate([experiment_scores, control_scores])
    cdf_experiment = np.searchsorted(experiment_scores, data_all, side='right') / en1
    cdf_control = np.searchsorted(control_scores, data_all, side='right') / en2

    # Compute the maximum distance between the two CDFs
    D = np.max(np.abs(cdf_experiment - cdf_control))

    # Calculate the KS statistic
    en = np.sqrt(en1 * en2 / (en1 + en2))
    p_value = ks((en + 0.12 + 0.11 / en) * D)

    return p_value


def ks(alam):
    """
    Compute the Kolmogorov-Smirnov probability given a lambda value.

    Parameters:
    - alam (float): Lambda value for the KS statistic.

    Returns:
    float: The probability associated with the KS statistic.
    """
    EPS1 = 1e-6  # Precision for the convergence of term's absolute value
    EPS2 = 1e-10  # Precision for the convergence of the series_sum's relative value
    a2 = -2.0 * alam ** 2  # Adjust lambda for exponential calculation
    fac = 2.0
    series_sum = 0.0
    previous_term = 0.0

    # Sum the series until convergence criteria are met
    for j in range(1, 101):
        term = fac * np.exp(a2 * j ** 2)  # Calculate term of the series
        series_sum += term  # Add to series_sum

        # Check for convergence
        if np.abs(term) <= EPS1 * previous_term or np.abs(term) <= EPS2 * series_sum:
            return series_sum

        fac = -fac  # Alternate the sign
        previous_term = np.abs(term)  # Update the term before flag

    # Return 1.0 if the series does not converge within 100 terms
    return 1.0


def compute_mw_python(experiment_ranks, control_ranks):
    """
    Compute the Mann-Whitney U test manually using rank sums to determine the statistical difference
    between two independent samples.

    Parameters:
    - experiment_ranks (list): Ranks of the experimental group.
    - control_ranks (list): Ranks of the control group.
    Returns:
    tuple: The Mann-Whitney U statistic and the corresponding p-value.
    """

    # Calculate the sum of ranks for each group
    R1 = np.sum(experiment_ranks)
    R2 = np.sum(control_ranks)

    # Number of observations in each group
    n1 = len(experiment_ranks)
    n2 = len(control_ranks)

    # Compute the Mann-Whitney U statistics for both groups
    U1 = R1 - n1 * (n1 + 1) / 2  # U statistic for the experimental group
    U2 = R2 - n2 * (n2 + 1) / 2  # U statistic for the control group

    # Use the smaller U statistic as the test statistic
    U = min(U1, U2)

    # Calculate the mean and standard deviation for the U distribution under H0 (null hypothesis)
    mean_U = n1 * n2 / 2
    std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

    # Calculate the Z-score associated with U statistic
    Z = (U - mean_U) / std_U

    # Calculate the two-tailed p-value from the Z-score
    p_value = 2 * norm.cdf(-np.abs(Z))  # Two-tailed test

    return U, p_value


def jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union


def hypergeometric_sf(x, M, N, n):
    """
    Calculate the survival function (complement of the CDF) for the hypergeometric distribution.

    Parameters:
    - x (int): The number of successful draws.
    - M (int): The total size of the population.
    - N (int): The number of success states in the population.
    - n (int): The number of draws.

    Returns:
    float: The survival function probability.
    """
    # Calculate the survival function using hypergeometric distribution
    probability = hypergeom.sf(x - 1, M, N, n)
    return probability