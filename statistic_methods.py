import numpy as np
from scipy.stats import ranksums, rankdata

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