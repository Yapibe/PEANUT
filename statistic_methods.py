import numpy as np
from scipy.stats import ranksums, rankdata, ttest_ind


class StatResults:
    def __init__(self, p_value=None, directionality=None, z_score=None, name=None):
        self.p_value = p_value
        self.directionality = directionality
        self.z_score = z_score
        self.name = name


def wilcoxon_rank_sums_test(experiment_scores, elements_scores, alternative='two-sided', **kwargs) -> StatResults:
    wilcoxon_rank_sums_test.name = 'Wilcoxon_rank_sums_test'
    p_vals = ranksums(experiment_scores, elements_scores, alternative=alternative).pvalue
    direction = np.mean(experiment_scores) > np.mean(elements_scores)
    return StatResults(p_value=p_vals, directionality=direction, name=wilcoxon_rank_sums_test.name)


from scipy.stats import ttest_ind
import numpy as np


def students_t_test(experiment_scores, elements_scores):
    """
    Perform a Student's t-test or Welch's t-test on two sets of data.

    :param experiment_scores: Scores from the experiment.
    :param elements_scores: Scores from the elements (population).
    :return: t-statistic, p-value, and directionality of the difference.
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
    p_vals_rank = rankdata(p_values, 'max') - 1
    p_vals_rank_ord = rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) / (p_vals_rank + 1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_vals_rank]
    return adj_p_vals
