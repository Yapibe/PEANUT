import numpy as np
from scipy.stats import ranksums


class StatRes:
    def __init__(self, p_value=None, directionality=None, z_score=None, name=None):
        self.p_value = p_value
        self.directionality = directionality
        self.z_score = z_score
        self.name = name


def wilcoxon_rank_sums_test(experiment_scores, elements_scores, alternative='two-sided', **kwargs) -> StatRes:
    wilcoxon_rank_sums_test.name = 'Wilcoxon_rank_sums_test'
    p_vals = ranksums(experiment_scores, elements_scores, alternative=alternative).pvalue
    direction = np.mean(experiment_scores) > np.mean(elements_scores)
    return StatRes(p_value=p_vals, directionality=direction, name=wilcoxon_rank_sums_test.name)
