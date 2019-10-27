import numpy as np
import scipy.special
import scipy.optimize


##############################################################################################################
##############################################################################################################
# This code has the functions used to plots that compare the optimal DP comp bounds from Murtagh and Vadhan '16
#  and the bounds for bounded range algos derived in the following  paper: https://arxiv.org/abs/1905.04273
# (where plots can be found).  The code to generate the plots that use the code here are in compare_comps.py.
##############################################################################################################
##############################################################################################################


# Version of DP optimal composition bound when each algorithm is (eps, delta = 0)-DP for the homogeneous case.
# We use the formula version from Jack Murtagh and Salil Vadhan's paper: https://arxiv.org/abs/1507.03113

def opt_compLHS(eps_g, eps, k):
    ell_start = int(np.ceil((eps_g + k * eps) / (2 * eps)))
    coeff = 1.0 / (1 + np.exp(eps)) ** k
    array_terms = []
    for i in range(ell_start, int(k + 1)):
        binom = scipy.special.binom(k, i)
        factor = np.exp(i * eps) - np.exp(eps_g) * np.exp((k - i) * eps)
        array_terms.append(binom * factor)
    return coeff * sum(array_terms)


# Use a root solver to find the eps_g that attains the delta bound.
def opt_compMV(k, eps, delta, eps_max=100):
    def f(eps_g):
        return opt_compLHS(eps_g, eps, k) - delta

    root = scipy.optimize.brentq(f, 0, eps_max)
    return root


# Now we use the formula for DP composition when each mechanism is
# epsilon-bounded range from this paper: https://arxiv.org/abs/1905.04273v2
def advanced_comp_br(k, eps, delta_thresh):
    if delta_thresh <= 0.0:
        return k * eps
    mean = k * eps ** 2 / 2.0
    sd = eps * np.sqrt(1.0 / 2.0 * k * np.log(1.0 / delta_thresh))
    mean_old = k * eps * (np.exp(eps) - 1.0) / (np.exp(eps) + 1.0)
    sd_old = eps * np.sqrt(2.0 * k * np.log(1.0 / delta_thresh))
    return min(mean + sd, mean_old + sd_old, k * eps)
